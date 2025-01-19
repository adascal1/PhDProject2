clc; clear all; close all; 
addpath(genpath('/mfip/mfip1/arielle/software/violin'))
%% (1) Distribution of Effect Size Across Metrics (i.e.: distribution of + or - effect size values across ROIs considered as hubs
% across the metrics) 
    hubPath = '/mfip/mfip1/arielle/PhDProject2/maps/MICs/hubLocations/';
    effectPath = '/mfip/mfip1/arielle/PhDProject2/maps/MICs/consistency/';

parcelNums = [100, 300, 600, 900]; 
group = 'MICs'

for iparcel = 1:length(parcelNums)

    % Loading the hub map ROIs
    Hub_lh = gifti(fullfile(hubPath, sprintf('schaefer%d_avgHubMap_betweenMetrics_thresh-0.100000.L.func.gii', parcelNums(iparcel)))).cdata;
    Hub_rh = gifti(fullfile(hubPath, sprintf('schaefer%d_avgHubMap_betweenMetrics_thresh-0.100000.R.func.gii', parcelNums(iparcel)))).cdata;

    hubs = [Hub_lh; Hub_rh];

    Effect_lh = gifti(fullfile(effectPath, sprintf('schaefer%d_avgCohenD_intraParcel_thresh-0.100000.L.func.gii', parcelNums(iparcel)))).cdata;
    Effect_rh = gifti(fullfile(effectPath, sprintf('schaefer%d_avgCohenD_intraParcel_thresh-0.100000.R.func.gii', parcelNums(iparcel)))).cdata;

    effect = [Effect_lh; Effect_rh];

    hubIdx = find(hubs);

    effectsAtHubs = unique(effect(hubIdx));

    T = table(repmat({sprintf('schaefer%d', parcelNums(iparcel))}, size(effectsAtHubs, 1), 1), effectsAtHubs, 'VariableNames', {'SchaeferParcellationNumber', 'UniqueEffectSizeValues'});
    writetable(T, sprintf('/mfip/mfip1/arielle/PhDProject2/results/%s/hubLocations/uniqueEffectSizesAtHubs.csv', group), 'WriteMode', 'append')

end
    
