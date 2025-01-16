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

% Creating new directories
if ~exist(sprintf('/mfip/mfip1/arielle/PhDProject2/code/%s/%s%d/richClub/', groupName, parcelName, parcelNum), 'dir')
    mkdir(sprintf('/mfip/mfip1/arielle/PhDProject2/code/%s/%s%d/richClub/', groupName, parcelName, parcelNum))
end

if ~exist(sprintf('/mfip/mfip1/arielle/PhDProject2/figures/%s/%s%d/richClub/', groupName, parcelName, parcelNum), 'dir')
    mkdir(sprintf('/mfip/mfip1/arielle/PhDProject2/figures/%s/%s%d/richClub/', groupName, parcelName, parcelNum))
end

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

    corrMatrix(:, :, iSub) = load(sprintf('/mfip/mfip1/arielle/PhDProject2/data/%s/%s%d_correlationMatrixes/sub-%s_ses-01_surf-fsLR-32k_parc-%s%d_desc-corrMatrix.mat', ...
        groupName, parcelName, parcelNum, sub_list{iSub}, parcelName, parcelNum)).corrMatrix;

end

avgCorrMatrix = mean(corrMatrix,3);

%% Applying Graph Metrics

% setting thresholds we want to do
 thresh = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]; 
for ithresh = 1:length(thresh)

    % applying weighted thresholding on the matrix
    threshMatrix = threshold_proportional(avgCorrMatrix, thresh(ithresh));

    % Computing Node Degree: 
    degree = sum((threshMatrix~=0)); %define degree of each node

    % computing the modualarity - the output 'modularity' is the community association vector.
    richClub = rich_club_wu(threshMatrix);

      
    % Reorder in a meaningfull maner 
    N = size(richClub, 2);
    idx = N - (0:N-1);
    richClub      = 1 - richClub(idx);

    % pulling out the rich club coefficient and the index it is located at (the index is the optimal node number)
    [richClubCoefficient, indexOfRichClubCoefficient] = max(richClub);


    % plotting the rich club curve (true value) against the null model rich club curves.
    % curve should be above all the null curves  
    figure; 
    hold on 
    p = plot(richClub, 'Color', [166/255,179/255,255/255], 'LineWidth', 3)
    xlabel('Degree Threshold (k)')
    ylabel('Rich Club Coefficient')
    title(sprintf('Rich Club Coefficient at Sparsity of %f', thresh(ithresh)))
    xlim([1, size(richClub, 2)])
    hold off
    
    % Add annotation with an arrow
    STR = sprintf('Rich Club found with %d nodes', indexOfRichClubCoefficient);
    annotation('textarrow',[0.4,0.4],[0.88,0.8],'String', STR,'TextColor','k','Fontsize', 16,'Color','w','LineStyle','-','FontWeight','bold');
    f = gcf;
    f.WindowState = 'Maximized';
    saveas(f, char(sprintf('/mfip/mfip1/arielle/PhDProject2/figures/%s/%s%d/richClub/richClubCurve_thresh-%f_group-%s.svg', ...
        groupName, parcelName, parcelNum, thresh(ithresh), subjType)))

    % Sorting the degrees and indices
    [degreesSort, nodeNumSort] = sort(degree, 'descend');

    % saving the degrees the size of our indexOfRichClubCoefficient
    idxOfNodesApartOfRichClub = nodeNumSort(1:indexOfRichClubCoefficient);    
    
    % plotting the rich club nodes on the cortex    
    toVisualize = zeros(size(parcelScheme));
    for idx = 1:length(idxOfNodesApartOfRichClub)
        
         idxBig = find(parcelScheme == idxOfNodesApartOfRichClub(idx));
         toVisualize(idxBig) = 1;

    end
   
   plot_hemispheres(toVisualize.*[mask_lh; mask_rh], {surf_lh, surf_rh}); colormap(slanCM('binary'))
   f = gcf;
   saveas(f, char(sprintf('/mfip/mfip1/arielle/PhDProject2/figures/%s/%s%d/richClub/richClubCortex_thresh-%f_group-%s.svg', ...
        groupName, parcelName, parcelNum, thresh(ithresh), subjType)))

    % Saving results - rich club coefficient
    T = table(richClub', 'VariableNames', {'richClub'});
    writetable(T, sprintf('/mfip/mfip1/arielle/PhDProject2/code/%s/%s%d/richClub/richClubCoefficientResults_thresh-%f_group-%s.csv', ...
        groupName, parcelName, parcelNum, thresh(ithresh), subjType));

    % saving results - nodes associated with rich club
    T = table(idxOfNodesApartOfRichClub', degreesSort(idxOfNodesApartOfRichClub)', 'VariableNames', {'idxOfNodesApartOfRichClub', 'degreesOfNodesApartOfRichClub'});
    writetable(T, sprintf('/mfip/mfip1/arielle/PhDProject2/code/%s/%s%d/richClub/richClubNodesResults_thresh-%f_group-%s.csv', ...
        groupName, parcelName, parcelNum, thresh(ithresh), subjType));
    
    clear richClub
    clear richClub_null
    
end

close all
