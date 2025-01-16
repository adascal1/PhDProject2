clc; clear all; close all;

%% Loading Data
path = '/mfip/mfip1/arielle/PhDProject2/';
addpath(genpath('/mfip/mfip1/arielle/software/matlab_GIfTI'));

% User input for analysis
parcellation = strsplit(input("What parcellation scheme and subject database do you want to do and metric? (i.e.: schaefer 100 MICs 0.1) \n" , "s"));
    parcelName  = parcellation{1};
    parcelNum   = str2double(parcellation{2});
    groupName   = parcellation{3};
    thresh      = str2double(parcellation{4});

% Loading parcellation
parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum)));

% Loading the cohenD results
DC = table2array(readtable(fullfile(path, 'results', groupName, 'degree', sprintf('degree_%s%d_cohendResults_thresh-%f.csv', parcelName, parcelNum, thresh))));
DC = DC(:, 1);
BC = table2array(readtable(fullfile(path, 'results', groupName, 'betweennessCentrality', sprintf('betweennessCentrality_%s%d_cohendResults_thresh-%f.csv', parcelName, parcelNum, thresh))));
BC = BC(:, 1);
EC = table2array(readtable(fullfile(path, 'results', groupName, 'eigenvectorCentrality', sprintf('eigenvectorCentrality_%s%d_cohendResults_thresh-%f.csv', parcelName, parcelNum, thresh))));
EC = EC(:, 1);

%% Finding where we have increases for all three metrics
% Finding just increases
DC_bin = DC > 0;
EB_bin = BC > 0;
EC_bin = EC > 0;

% Computing the Intersection over union for above to calculate the similarity 
% Compute the intersection (common active regions in all three maps)
intersection_map_above = DC_bin & EB_bin & EC_bin;
union_map_above        = DC_bin | EB_bin | EC_bin;

% computing the intersection over union 
iou_above = sum(intersection_map_above)/sum(union_map_above);

%% Finding where we have decreases for all three metrics
% Finding just increases
DC_bin = DC < 0;
EB_bin = BC < 0;
EC_bin = EC < 0;

% Computing the Intersection over union for above to calculate the similarity 
% Compute the intersection (common active regions in all three maps)
intersection_map_below = DC_bin & EB_bin & EC_bin;
union_map_below        = DC_bin | EB_bin | EC_bin;

% computing the intersection over union 
iou_below = sum(intersection_map_below)/sum(union_map_below);

%% Writing A table
T = table({sprintf('%s%d', parcelName, parcelNum)}, iou_above, iou_below, 'VariableNames', {'Parcellation', 'iouAbove', 'iouBelow'});
writetable(T, fullfile(path, 'results', 'MICs', 'consistency', 'intersectionOverUnion_intraParcellation.csv'), 'WriteMode', 'append')
