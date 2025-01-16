clc; clear all; close all;

%% Loading Data
path = '/mfip/mfip1/arielle/PhDProject2/';
addpath(genpath('/mfip/mfip1/arielle/software/matlab_GIfTI'));

% User input for analysis
parcellation = strsplit(input("What parcellation schemes? (i.e.: 100 300) \n" , "s"));
    parcelName   = 'schaefer'; %parcellation{1};
    parcelNum1   = str2double(parcellation{1});
    parcelNum2   = str2double(parcellation{2});
    groupName    = 'MICs'; %parcellation{3};
    thresh       = 0.1; %str2double(parcellation{4});

% Loading parcellation
parcelScheme1 = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum1)));
parcelScheme2 = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum2)));

% Loading the cohenD results from the maps to study at a vertex level- parcel1 (LOWER NUMBER PARCELATION)
DC_parcel1 = load_cohendMaps(groupName, 'degree', parcelName, parcelNum1, thresh);
BC_parcel1 = load_cohendMaps(groupName, 'betweennessCentrality', parcelName, parcelNum1, thresh);
EC_parcel1 = load_cohendMaps(groupName, 'eigenvectorCentrality', parcelName, parcelNum1, thresh);

% Loading the cohenD results from the maps to study at a vertex level- parcel2 (HIGHER NUMBER PARCELATION)
DC_parcel2 = load_cohendMaps(groupName, 'degree', parcelName, parcelNum2, thresh);
BC_parcel2 = load_cohendMaps(groupName, 'betweennessCentrality', parcelName, parcelNum2, thresh);
EC_parcel2 = load_cohendMaps(groupName, 'eigenvectorCentrality', parcelName, parcelNum2, thresh);

% Goal : get a similarity value for each schefer100 regions when compared to the larger parcellations

%% IOU
metrics = {'degree', 'betweennessCentrality', 'eigenvectorCentrality'};
for iMetric = 1:length(metrics)
    if strcmp(metrics{iMetric}, 'degree')
        parcel1 = DC_parcel1;
        parcel2 = DC_parcel2;
    elseif strcmp(metrics{iMetric}, 'betweennessCentrality')
        parcel1 = BC_parcel1;
        parcel2 = BC_parcel2;
    elseif strcmp(metrics{iMetric}, 'eigenvectorCentrality')
        parcel1 = EC_parcel1;
        parcel2 = EC_parcel2;
    end

    for iParcel1 = 1:(length(unique(parcelScheme1))-1)
            
        % Fining the indices of the iParcel1 number association with the LARGER parcellation scheme        
        idxParcel1 = find(parcelScheme1 == iParcel1);
        parcel1Val = unique(parcel1(idxParcel1));
        if parcel1Val > 0
            parcel1_bin = parcel1 > 0;
            parcel2_bin = parcel2 > 0;
            % If we have a positive or negative effect size
            posOrneg = 'positive';
        elseif parcel1Val < 0 
            parcel1_bin = parcel1 < 0;
            parcel2_bin = parcel2 < 0;
            posOrneg = 'negative';
        end
        
        %% Pulling out the values at both parcels
        parcel1_bin_vals= parcel1_bin(idxParcel1);
        parcel2_bin_vals= parcel2_bin(idxParcel1);
    
        % computing IOU for that parcel:
        intersection_map = parcel1_bin_vals & parcel2_bin_vals;
        union_map        = parcel1_bin_vals | parcel1_bin_vals;
    
        iou              = sum(intersection_map)/sum(union_map);
    
        T = table({parcelName}, parcelNum1, parcelNum2, {sprintf('Region%d', iParcel1)}, {metrics{iMetric}}, iou, {posOrneg}, ...
            'VariableNames', {'ParcelationName', 'ParcelNumber1', 'ParcelNumber2', 'ParcelNumber1Region', 'metric', 'intersectionOverUnion', 'PositiveOrNegativeEffectSize'});

        writetable(T, fullfile(path, 'results', 'MICs', 'consistency', 'intersectionOverUnion_interParcellation.csv'), 'WriteMode', 'append');
    end

end





















function [metricOut] = load_cohendMaps(groupName, metric, parcelName, parcelNum, thresh)

path = '/mfip/mfip1/arielle/PhDProject2/';

tmp1 = gifti(fullfile(path, 'maps', groupName, metric, sprintf('%s_%s%d_thresh-%f_cohend.L.func.gii', metric, parcelName, parcelNum, thresh))).cdata;
tmp2 = gifti(fullfile(path, 'maps', groupName, metric, sprintf('%s_%s%d_thresh-%f_cohend.R.func.gii', metric, parcelName, parcelNum, thresh))).cdata;
metricOut = [tmp1; tmp2];

end
