clc; clear all; close all; 

%% Loading paths
addpath(genpath('/mfip/mfip1/arielle/software/BrainConnectivityToolbox'));

% Getting user input for parcellation scheme to use
parcellation = strsplit(input("What parcellation scheme and subject database do you want to do? (i.e.: schaefer 100 MICs HC) \n" , "s"));
    parcelName = parcellation{1};
    parcelNum  = str2double(parcellation{2});
    groupName = parcellation{3};
    subjType = parcellation{4};

% Creating new directories
if ~exist(sprintf('/mfip/mfip1/arielle/PhDProject2/data/%s/clusteringCoefficient/%s%d/', groupName, parcelName, parcelNum), 'dir')
    mkdir(sprintf('/mfip/mfip1/arielle/PhDProject2/data/%s/clusteringCoefficient/%s%d/', groupName, parcelName, parcelNum))
end

%% Loading Subject Data
if strcmp(subjType, 'PX')
    info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_list_%s_%s_final.xlsx', groupName, subjType));
    sub_list = info.Subj_ID;
elseif strcmp(subjType, 'HC')
    info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_list_%s_%s_final.xlsx', groupName, subjType));
    sub_list = info.Subj_ID;
end

% setting thresholds we want to do
thresh = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5];
for iSub = 1:size(info, 1)

     corrMatrix = load(sprintf('/mfip/mfip1/arielle/PhDProject2/data/%s/correlationMatrixes/%s%d_correlationMatrixes/sub-%s_ses-01_surf-fsLR-32k_parc-%s%d_desc-corrMatrix.mat', ...
        groupName, parcelName, parcelNum, sub_list{iSub}, parcelName, parcelNum)).corrMatrix;


     avgCorrMatrix = mean(corrMatrix,3);


    %% Applying Graph Metrics
    
    for ithresh = 1:length(thresh)
    
        % applying weighted thresholding on the matrix
        threshMatrix = threshold_proportional(avgCorrMatrix, thresh(ithresh));
    
        % computing the eigenvector centrality
        clusteringCoefficient = clustering_coef_wu(threshMatrix)';
        
    
        % Saving results
        T = table(clusteringCoefficient', 'VariableNames', {'clusteringCoefficient'});
        writetable(T, sprintf('/mfip/mfip1/arielle/PhDProject2/data/%s/clusteringCoefficient/%s%d/clusteringCoefficientResults_sub-%s_thresh-%f.csv', ...
            groupName, parcelName, parcelNum, sub_list{iSub}, thresh(ithresh)));
            
    end
end