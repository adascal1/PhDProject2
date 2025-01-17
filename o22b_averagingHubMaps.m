clc; clear all; close all; 

path = '/mfip/mfip1/arielle/PhDProject2/';
addpath(genpath('/mfip/mfip1/arielle/software/matlab_GIfTI'));

    parcelName  = 'schaefer'; %parcellation{1};
    groupName   = 'MICs'; %parcellation{3};
    thresh      = 0.1; %str2double(parcellation{4});

% Loading parcellation
%parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum)));

parcelNums = {'100', '300', '600', '900'};

%% Avereraging across the metrics (to determine the strength of the hub)

for iparcelNum = 1:length(parcelNums)

    % Degree Centrality
    lh_DC = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s%s_metric-degreeCentrality_thresh-%f_subjType-HC.L.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;
    rh_DC = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s%s_metric-degreeCentrality_thresh-%f_subjType-HC.R.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;

    % Betweenness Centrality
    lh_BC = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s%s_metric-betweennessCentrality_thresh-%f_subjType-HC.L.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;
    rh_BC = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s%s_metric-betweennessCentrality_thresh-%f_subjType-HC.R.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;
    
    % Participation Coefficient
    lh_PC = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s%s_metric-participationCoefficient_thresh-%f_subjType-HC.L.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;
    rh_PC = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s%s_metric-participationCoefficient_thresh-%f_subjType-HC.R.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;

    lh = [lh_DC, lh_BC, lh_PC];
    rh = [rh_DC, rh_BC, rh_PC];
    
    mean_lh = mean(lh, 2);
    mean_rh = mean(rh, 2);

    savegifti(mean_lh, mean_rh, groupName, parcelName, parcelNums{iparcelNum}, thresh, 'betweenMetrics');

    
end

%% Averaging across the parcels (to determine the consistentcy of the hub)
metrics = {'degreeCentrality', 'betweennessCentrality', 'participationCoefficient'};

for imetric = 1:length(metrics)

    lh_100 = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s100_metric-%s_thresh-%f_subjType-HC.L.func.gii', parcelName, metrics{imetric}, thresh))).cdata;
    rh_100 = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s100_metric-%s_thresh-%f_subjType-HC.R.func.gii', parcelName, metrics{imetric}, thresh))).cdata;

    lh_300 = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s300_metric-%s_thresh-%f_subjType-HC.L.func.gii', parcelName, metrics{imetric}, thresh))).cdata;
    rh_300 = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s300_metric-%s_thresh-%f_subjType-HC.R.func.gii', parcelName, metrics{imetric}, thresh))).cdata;

    lh_600 = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s600_metric-%s_thresh-%f_subjType-HC.L.func.gii', parcelName, metrics{imetric}, thresh))).cdata;
    rh_600 = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s600_metric-%s_thresh-%f_subjType-HC.R.func.gii', parcelName, metrics{imetric}, thresh))).cdata;

    lh_900 = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s900_metric-%s_thresh-%f_subjType-HC.L.func.gii', parcelName, metrics{imetric}, thresh))).cdata;
    rh_900 = gifti(fullfile(path, 'maps', groupName, 'hubLocations', sprintf('hubLocations_%s900_metric-%s_thresh-%f_subjType-HC.R.func.gii', parcelName, metrics{imetric}, thresh))).cdata;

    lh = [lh_100, lh_300, lh_600, lh_900];
    rh = [rh_100, rh_300, rh_600, rh_900];
    
    mean_lh = mean(lh, 2);
    mean_rh = mean(rh, 2);

    savegifti(mean_lh, mean_rh, groupName, parcelName, parcelNums{iparcelNum}, thresh, metrics{imetric});

    
end







%% Extra Functions
function savegifti(mapToSaveLeft, mapToSaveRight, groupName, parcelName, parcelNum, thresh, whichAvg)

if strcmp(whichAvg, 'betweenMetrics')
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/hubLocations/%s%s_avgHubMap_%s_thresh-%f.L.func.gii', groupName, parcelName, parcelNum, whichAvg, thresh))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/hubLocations/%s%s_avgHubMap_%s_thresh-%f.L.func.gii', groupName, parcelName, parcelNum, whichAvg, thresh))
else 
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/hubLocations/%s_avgHubMap_%s_thresh-%f.L.func.gii', groupName, parcelName, whichAvg, thresh))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/hubLocations/%s_avgHubMap_%s_thresh-%f.R.func.gii', groupName, parcelName, whichAvg, thresh))
end

end
