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
toVisualize = zeros(size(parcelScheme1));

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

    % computing the RMSE between each parcel pairing
    for iParcel1 = 1:(length(unique(parcelScheme1))-1)
            
        idx = find(parcelScheme1 == iParcel1);
        parcel1_idx = parcel1(idx); 
        parcel2_idx = parcel2(idx);

        diff = (parcel1_idx - parcel2_idx).^2; % difference between each vertex squared
        SSE = sum(diff); % sum of the differences between the vertices
        RMSE = sqrt(SSE/length(idx)); % number of vertices

                  
        T = table({parcelName}, parcelNum1, parcelNum2, {sprintf('Region%d', iParcel1)}, {metrics{iMetric}}, SSE, RMSE, ...
            'VariableNames', {'ParcelationName', 'ParcelNumber1', 'ParcelNumber2', 'ParcelNumber1Region', 'metric', 'SSE', 'RMSE'});

        writetable(T, fullfile(path, 'results', 'MICs', 'consistency', 'SSE_RMSE_interParcellation.csv'), 'WriteMode', 'append');

        toVisualize(idx) = RMSE; % placing the RMSE at these index values
        % splitting the left and right hemispheres for plotting
        lh = toVisualize(1:size(toVisualize, 1)/2);
        rh = toVisualize((size(toVisualize, 1)/2)+1:end);

        % saving the gifti files
        savegifti(lh, rh, groupName, metrics{iMetric}, parcelName, parcelNum1, parcelNum2, thresh)



    end

end





%% Extra Functions
function [metricOut] = load_cohendMaps(groupName, metric, parcelName, parcelNum, thresh)

path = '/mfip/mfip1/arielle/PhDProject2/';

tmp1 = gifti(fullfile(path, 'maps', groupName, metric, sprintf('%s_%s%d_thresh-%f_cohend.L.func.gii', metric, parcelName, parcelNum, thresh))).cdata;
tmp2 = gifti(fullfile(path, 'maps', groupName, metric, sprintf('%s_%s%d_thresh-%f_cohend.R.func.gii', metric, parcelName, parcelNum, thresh))).cdata;
metricOut = [tmp1; tmp2];

end

function savegifti(mapToSaveLeft, mapToSaveRight, groupName, metric, parcelName, parcelNum1, parcelNum2, thresh)
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/%s_%s%dvs%d_thresh-%f_RMSE.L.func.gii', groupName, metric, parcelName, parcelNum1, parcelNum2, thresh))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/%s_%s%dvs%d_thresh-%f_RMSE.R.func.gii', groupName, metric, parcelName, parcelNum1, parcelNum2, thresh))
end

