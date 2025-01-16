clc; clear all; close all; 

path = '/mfip/mfip1/arielle/PhDProject2/';
addpath(genpath('/mfip/mfip1/arielle/software/matlab_GIfTI'));

% User input for analysis
%parcellation = strsplit(input("What parcellation scheme and subject database do you want to do and metric? (i.e.: schaefer 100 MICs 0.1) \n" , "s"));
    parcelName  = 'schaefer'; %parcellation{1};
    % parcelNum   = 100; %str2double(parcellation{2});
    groupName   = 'MICs'; %parcellation{3};
    thresh      = 0.1; %str2double(parcellation{4});

% Loading parcellation
%parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum)));

metrics = {'degree', 'betweennessCentrality', 'eigenvectorCentrality', 'participationCoefficient'};

%% Loading the Maps

for imetric = 1:length(metrics)

    lh_100 = gifti(fullfile(path, 'maps', groupName, metrics{imetric}, sprintf('%s_%s100_thresh-%f_cohend.L.func.gii', metrics{imetric}, parcelName, thresh))).cdata;
    rh_100 = gifti(fullfile(path, 'maps', groupName, metrics{imetric}, sprintf('%s_%s100_thresh-%f_cohend.R.func.gii', metrics{imetric}, parcelName, thresh))).cdata;

    lh_300 = gifti(fullfile(path, 'maps', groupName, metrics{imetric}, sprintf('%s_%s300_thresh-%f_cohend.L.func.gii', metrics{imetric}, parcelName, thresh))).cdata;
    rh_300 = gifti(fullfile(path, 'maps', groupName, metrics{imetric}, sprintf('%s_%s300_thresh-%f_cohend.R.func.gii', metrics{imetric}, parcelName, thresh))).cdata;

    lh_600 = gifti(fullfile(path, 'maps', groupName, metrics{imetric}, sprintf('%s_%s600_thresh-%f_cohend.L.func.gii', metrics{imetric}, parcelName, thresh))).cdata;
    rh_600 = gifti(fullfile(path, 'maps', groupName, metrics{imetric}, sprintf('%s_%s600_thresh-%f_cohend.R.func.gii', metrics{imetric}, parcelName, thresh))).cdata;

    lh_900 = gifti(fullfile(path, 'maps', groupName, metrics{imetric}, sprintf('%s_%s900_thresh-%f_cohend.L.func.gii', metrics{imetric}, parcelName, thresh))).cdata;
    rh_900 = gifti(fullfile(path, 'maps', groupName, metrics{imetric}, sprintf('%s_%s900_thresh-%f_cohend.R.func.gii', metrics{imetric}, parcelName, thresh))).cdata;

    lh = [lh_100, lh_300, lh_600, lh_900];
    rh = [rh_100, rh_300, rh_600, rh_900];
    
    mean_lh = mean(lh, 2);
    mean_rh = mean(rh, 2);

    savegifti(mean_lh, mean_rh, groupName, metrics{imetric}, thresh);

    
end




%% Extra Functions
function savegifti(mapToSaveLeft, mapToSaveRight, groupName, metric, thresh)
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/%s_avgCohenD_interParcel_thresh-%f.L.func.gii', groupName, metric, thresh))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/%s_avgCohenD_interParcel_thresh-%f.R.func.gii', groupName, metric, thresh))
end