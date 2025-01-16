clc; clear all; close all;
addpath(genpath('/mfip/mfip1/arielle/software/BrainConnectivityToolbox'));
addpath(genpath('/mfip/mfip1/arielle/software/fdr_bh'));

%% For Visualization
addpath(genpath('/mfip/mfip1/arielle/software/BrainSpace'));
addpath(genpath('/mfip/mfip1/arielle/software/matlab_GIfTI'));
addpath(genpath('/mfip/mfip1/arielle/software/slanCM'));

[surf_lh, surf_rh] = load_conte69();
[mask_lh, mask_rh] = load_mask();

%% Getting user input for parcellation scheme to use
parcellation = strsplit(input("What parcellation scheme and subject database do you want to do? (i.e.: schaefer 100 MICs HC) \n" , "s"));
    parcelName = parcellation{1};
    parcelNum  = str2double(parcellation{2});
    groupName = parcellation{3};
    subjType = parcellation{4};

%% Getting user input for parcellation scheme to use
parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum)));
toVisualize = zeros(size(parcelScheme));    

%% Create Average Matrix
if strcmp(subjType, 'PX')
    info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_list_%s_final.xlsx', groupName));
    sub_list = info.Subj_ID;
elseif strcmp(subjType, 'HC')
    info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_list_%s_HC_final.xlsx', groupName));
    sub_list = info.Subj_ID;
end

for iSub = 1:size(info, 1)

  corrMatrix(:, :, iSub) = load(sprintf('/mfip/mfip1/arielle/PhDProject2/data/%s/correlationMatrixes/%s%d_correlationMatrixes/sub-%s_ses-01_surf-fsLR-32k_parc-%s%d_desc-corrMatrix.mat', ...
        groupName, parcelName, parcelNum, sub_list{iSub}, parcelName, parcelNum)).corrMatrix;

end

avgCorrMatrix = mean(corrMatrix,3);


%% Applying Graph Metrics

% setting thresholds we want to do
%thresh = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5];
thresh = 0.1;

for ithresh = 1:length(thresh)

    % applying weighted thresholding on the matrix
    threshMatrix = threshold_proportional(avgCorrMatrix, thresh(ithresh));

    % computing the eigenvector centrality
    eigenvectorCentrality = eigenvector_centrality_und(threshMatrix);
    % [m, I] = sort(eigenvectorCentrality,'descend'); % m = ordered eigenvector centrality values; I = ordered node index values
    % tmp = strsplit(sprintf('%d ', I));
    % 
    % % plotting the eigenvector centrality values, with a line at the top 10% of nodes
    % figure; 
    % bar(m(1:50),'FaceAlpha', 0.5,'FaceColor',[0.659    0.573    0.878]); 
    % xticks(1:1:length(m)); 
    % xticklabels(tmp(1:1:end-1))
    % thresNode = quantile(m, 0.9); % the top 10% nodes
    % yline(thresNode, '--', 'Color', [0.98 0.706, 0.733], 'LineWidth', 3); 
    % % xline( max(find(m>thresNode))+0.5,'--', 'Color', [0.98 0.706, 0.733], 'LineWidth', 3);
    % xlabel('Node Number')
    % ylabel('Eigenvector Centrality')
    % xlim([0.5 50.5])
    % title(sprintf('Sparsity of %f', thresh(ithresh)))
    % f = gcf;
    % f.WindowState = 'maximized';
    % saveas(f, char(sprintf('/mfip/mfip1/arielle/PhDProject2/figures/%s/%s%d/eigenvectorCentrality/eigenvectorCentralityBarplot_thresh-%f_group-%s.svg', ...
    %     groupName, parcelName, parcelNum, thresh(ithresh), subjType)));

    % if i want to put in a patch to show the roi ...
    % ytmp = ylim;
    % xtmp = xlim;
    % patch([xtmp(1) xtmp(2)/2 xtmp(2)/2 xtmp(1)],[0 0 ytmp(2) ytmp(2)],'red','FaceAlpha',0.2) 

    % Saving results
    T = table(eigenvectorCentrality, 'VariableNames', {'EigenVectorCentrality'});
        writetable(T, sprintf('/mfip/mfip1/arielle/PhDProject2/results/%s/eigenvectorCentrality/%s%d_averageEigenvectorCentralityResults_thresh-%f_group-%s.csv', ...
        groupName, parcelName, parcelNum, thresh(ithresh), subjType));
    
     
       
    % Visualization using brainspace - only visualizing nodes with significant p values
    toVisualize = zeros(size(parcelScheme));
    
    for idx = 1:size(avgCorrMatrix, 1)
         idxBig = find(parcelScheme == idx);
         toVisualize(idxBig) = eigenvectorCentrality(idx);
    end
   
    % Normalizing Results
    toVisualize = toVisualize/max(toVisualize);
    
    % splitting the left and right hemispheres for plotting
    lh = toVisualize(1:size(toVisualize, 1)/2);
    rh = toVisualize((size(toVisualize, 1)/2)+1:end);
    
    % saving the gifti files
    savegifti(lh, rh, groupName, 'eigenvectorCentrality', parcelName, parcelNum, thresh, subjType)
   
  
end

%% Extra functions
function savegifti(mapToSaveLeft, mapToSaveRight, groupName, metric, parcelName, parcelNum, thresh, subjType)
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/%s/%s_%s%d_thresh-%f_subjType-%s_averageMap.L.func.gii', groupName, metric, metric, parcelName, parcelNum, thresh, subjType))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/%s/%s_%s%d_thresh-%f_subjType-%s_averageMap.R.func.gii', groupName, metric, metric, parcelName, parcelNum, thresh, subjType))
end