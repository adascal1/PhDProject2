clc; 
clear all; 
close all;
%% This script identifies where hubs are located according to the following thresholds: 
% Participation Coefficient: Having a value > 0.6
% Betweenness Centrality : Taking the top 20% of nodes with the highest BC values, after removing the nodes with zero
% Degree Centrality: taking the top 20% of nodes with the highest DC value, and the bottom 20% of nodes with the lowest DC value, after removing the nodes with zero.

dataPath = '/mfip/mfip1/arielle/PhDProject2/results/MICs/';
threshold = 0.1;

%% Looping through parcellations: 

parcelNum = [100, 300, 600, 900];

for iParcel = 1:length(parcelNum)

    % Loading Parcel Scheme
    parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/schaefer-%d_conte69.csv', parcelNum(iParcel))));


    % participation coefficient
    PC = table2array(readtable(fullfile(dataPath, 'participationCoefficient', sprintf('schaefer%d_averageParticipationCoefficientResults_thresh-%f_group-HC.csv', parcelNum(iParcel), threshold))));
    hubs_PC_bin = PC > 0.6;
    hubs_PC = PC.*hubs_PC_bin;

    % betweenness centrality
    BC = table2array(readtable(fullfile(dataPath, 'betweennessCentrality', sprintf('schaefer%d_averageBetweennessCentralityResults_thresh-%f_group-HC.csv', parcelNum(iParcel), threshold))));
    thresh = quantile(BC(BC ~=0), 0.8); % finding the top 20% of nodes
    hubs_BC_bin = BC > thresh; 
    hubs_BC = BC.*hubs_BC_bin;
    
    % degree centrality
    DC = table2array(readtable(fullfile(dataPath, 'degree', sprintf('schaefer%d_averageDegreeResults_thresh-%f_group-HC.csv', parcelNum(iParcel), threshold))));
    threshs  = quantile(DC(DC ~=0 ), [0.2, 0.8]); % finding the bottom 20% of nodes and the top 20% of nodes
    thresh_low = threshs(1); thresh_high = threshs(2);
    hubs_DC_bin = DC < thresh_low | DC > thresh_high;
    hubs_DC = DC.*hubs_DC_bin;

    % Creating a table with all of this information
    T = table(hubs_PC, hubs_BC, hubs_DC, 'VariableNames', {'ParticipationCoefficientHubs', 'BetweennessCentralityHubs', 'DegreeCentralityHubs'});
    writetable(T, sprintf('/mfip/mfip1/arielle/PhDProject2/results/MICs/hubLocations/thresholdedHubMaps_schaefer%d_thresh-%f_group-HC.csv', ...
            parcelNum(iParcel), threshold));


    % Visualization using brainspace 
    metrics = {'betweennessCentrality', 'participationCoefficient', 'degreeCentrality'};

    for iMetric = 1:length(metrics)

        toVisualize = zeros(size(parcelScheme));
        toVisualize_bin = zeros(size(parcelScheme));


        if strcmp(metrics(iMetric), 'betweennessCentrality')
            data_bin = hubs_BC_bin;
            data     = hubs_BC;

        elseif strcmp(metrics(iMetric), 'participationCoefficient')
            data_bin = hubs_PC_bin;
            data     = hubs_PC;

        elseif strcmp(metrics(iMetric), 'degreeCentrality')
            data_bin = hubs_DC_bin;
            data     = hubs_DC;

        end

        for idx = 1:length(data_bin)
            idxBig = find(parcelScheme == idx);
            toVisualize_bin(idxBig) = data_bin(idx); % binarized maps
            toVisualize(idxBig) = data(idx); % value maps
            
        end
       
        % Normalization of Metrics
        toVisualize = toVisualize/max(toVisualize);

        % splitting the left and right hemispheres for plotting - binary
        lh_bin = toVisualize_bin(1:size(toVisualize_bin, 1)/2);
        rh_bin = toVisualize_bin((size(toVisualize_bin, 1)/2)+1:end);
        
        % saving the gifti files - binary
     %   savegifti(lh_bin, rh_bin, parcelNum(iParcel), sprintf('%sBin', metrics{iMetric}), threshold)
        
        % splitting the left and right hemispheres for plotting - values
        lh = toVisualize(1:size(toVisualize, 1)/2);
        rh = toVisualize((size(toVisualize, 1)/2)+1:end);
        
        % saving the gifti files
        savegifti(lh, rh, parcelNum(iParcel), sprintf('%s', metrics{iMetric}), threshold)


    end
    
      
end

%% Extra functions
function savegifti(mapToSaveLeft, mapToSaveRight, parcelNum, metric, thresh)
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/MICs/hubLocations/hubLocations_schaefer%d_metric-%s_thresh-%f_subjType-HC.L.func.gii', parcelNum, metric, thresh))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/MICs/hubLocations/hubLocations_schaefer%d_metric-%s_thresh-%f_subjType-HC.R.func.gii', parcelNum, metric, thresh))
end