clc; clear all; close all;

%% Loading paths
addpath(genpath('/mfip/mfip1/arielle/software/BrainConnectivityToolbox'));
addpath(genpath('/mfip/mfip1/arielle/software/BrainSpace'));
addpath(genpath('/mfip/mfip1/arielle/software/matlab_GIfTI'));
addpath(genpath('/mfip/mfip1/arielle/software/fdr_bh'));
addpath(genpath('/mfip/mfip1/arielle/software/slanCM'));

%% Loading visualization tools
[surf_lh, surf_rh] = load_conte69();
[mask_lh, mask_rh] = load_mask();
% Getting user input for parcellation scheme to use
parcellation = strsplit(input("What parcellation scheme and subject database do you want to do? (i.e.: schaefer 100 MICs HC) \n" , "s"));
    parcelName = parcellation{1};
    parcelNum  = str2double(parcellation{2});
    groupName = parcellation{3};
    subjType = parcellation{4};

%% Getting user input for parcellation scheme to use
parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum)));
toVisualize = zeros(size(parcelScheme));    

% Creating new directories
if ~exist(sprintf('/mfip/mfip1/arielle/PhDProject2/code/%s/%s%d/clusteringCoefficient/', groupName, parcelName, parcelNum), 'dir')
    mkdir(sprintf('/mfip/mfip1/arielle/PhDProject2/code/%s/%s%d/clusteringCoefficient/', groupName, parcelName, parcelNum))
end

if ~exist(sprintf('/mfip/mfip1/arielle/PhDProject2/figures/%s/%s%d/clusteringCoefficient/', groupName, parcelName, parcelNum), 'dir')
    mkdir(sprintf('/mfip/mfip1/arielle/PhDProject2/figures/%s/%s%d/clusteringCoefficient/', groupName, parcelName, parcelNum))
end


%% Create Average Matrix
if strcmp(subjType, 'PX')
    info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_list_%s_final.xlsx', groupName));
    sub_list = info.Subj_ID;
elseif strcmp(subjType, 'HC')
    info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_list_%s_HC_final.xlsx', groupName));
    sub_list = info.Subj_ID;
end

for iSub = 1:size(info, 1)

    corrMatrix(:, :, iSub) = load(sprintf('/mfip/mfip1/arielle/PhDProject2/data/%s/%s%d_correlationMatrixes/sub-%s_ses-01_surf-fsLR-32k_parc-%s%d_desc-corrMatrix.mat', ...
        groupName, parcelName, parcelNum, sub_list{iSub}, parcelName, parcelNum)).corrMatrix;

end

avgCorrMatrix = mean(corrMatrix,3);



%% Calculate Clustering Coefficient
% setting thresholds we want to do
thresh = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5] %, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1];
% thresh = [0.45, 0.5] %, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1];

for ithresh = 1:length(thresh)

    % Thresholding the matrix
    threshMatrix = threshold_proportional(avgCorrMatrix, thresh(ithresh));

    % Calculating clustering coefficient (sum of the binary matrix, number of connections found at that node)
    clusteringCoefficient = clustering_coef_wu(threshMatrix)';

    % Saving results
    T = table(clusteringCoefficient', 'VariableNames', {'ClusteringCoefficient'});
    writetable(T, sprintf('/mfip/mfip1/arielle/PhDProject2/code/%s/%s%d/clusteringCoefficient/clusteringCoefficientResults_thresh-%f_group-%s.csv', ...
        groupName, parcelName, parcelNum, thresh(ithresh), subjType));

    % Visualization using brainspace - only visualizing nodes with significant p values
    toVisualize = zeros(size(parcelScheme));
    for idx = 1:size(avgCorrMatrix, 1)
            idxBig = find(parcelScheme == idx);
            toVisualize(idxBig) = clusteringCoefficient(idx);
    end
    
    plot_hemispheres(toVisualize.*[mask_lh; mask_rh], {surf_lh, surf_rh}); colormap(slanCM('YlGnBu'))
    f = gcf;
    saveas(f, char(sprintf('/mfip/mfip1/arielle/PhDProject2/figures/%s/%s%d/clusteringCoefficient/clusteringCoefficient_thresh-%f_group-%s.svg', ...
        groupName, parcelName, parcelNum, thresh(ithresh), subjType)))

end

close all