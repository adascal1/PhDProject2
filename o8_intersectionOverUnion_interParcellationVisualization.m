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

data = readtable(fullfile(path, 'results', 'MICs', 'consistency', 'intersectionOverUnion_interParcellation.csv'));
parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum1)));

%% Getting information from the table to plot: 
% List of metrics
metrics = {'degree', 'betweennessCentrality', 'eigenvectorCentrality'};

% Data from just the two parcelations of interest
% subsectData = data(find(data.ParcelNumber1 == parcelNum1 & data.ParcelNumber2 == parcelNum2), :);

for iMetric = 1:length(metrics)

    % Data from just the two parcelations of interest and the metric of interest
    metricData = data(find(data.ParcelNumber1 == parcelNum1 & data.ParcelNumber2 == parcelNum2 & strcmp(data.metric, metrics{iMetric})), :);
    
    % creating gifti files
    toVisualize = zeros(size(parcelScheme));
    for idx = 1:size(metricData, 1)
    
        % If we had a negative effect size, setting the iou to negative to represent that
        if strcmp(metricData.PositiveOrNegativeEffectSize(idx), 'negative')
            metricData.intersectionOverUnion(idx) =  -metricData.intersectionOverUnion(idx);
        end

        idxBig = find(parcelScheme == idx); % finding the indexes in the parcel scheme associated with the parcel number 
        toVisualize(idxBig) = metricData.intersectionOverUnion(idx); % placing the cohend at these index values
    
    end

    % splitting the left and right hemispheres for plotting
    lh = toVisualize(1:size(toVisualize, 1)/2);
    rh = toVisualize((size(toVisualize, 1)/2)+1:end);

    % saving the gifti files
    savegifti(lh, rh, groupName, metrics{iMetric}, parcelName, parcelNum1, parcelNum2, thresh)


end



%% Extra Functions
function savegifti(mapToSaveLeft, mapToSaveRight, groupName, metric, parcelName, parcelNum1, parcelNum2, thresh)
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/%s_%s%dvs%d_thresh-%f_IOU.L.func.gii', groupName, metric, parcelName, parcelNum1, parcelNum2, thresh))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/%s_%s%dvs%d_thresh-%f_IOU.R.func.gii', groupName, metric, parcelName, parcelNum1, parcelNum2, thresh))
end