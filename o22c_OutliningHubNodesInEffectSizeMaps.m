clc; clear all; close all; 


%% Loading Schaefer 100 Hub Map
thresh = 0.1;
group = 'MICs';
parcelNum = [100, 300, 600, 900];

for iParcelNum = 1:length(parcelNum)

    % Loading the parcelScheme
    parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/schaefer-%d_conte69.csv', parcelNum(iParcelNum))));

    % Loading the Hub Maps averaged across metrics within the same parcellation
    hubMap_lh = gifti(sprintf('/mfip/mfip1/arielle/PhDProject2/maps/MICs/hubLocations/schaefer%d_avgHubMap_betweenMetrics_thresh-%f.L.func.gii', parcelNum(iParcelNum), thresh)).cdata;
    hubMap_rh = gifti(sprintf('/mfip/mfip1/arielle/PhDProject2/maps/MICs/hubLocations/schaefer%d_avgHubMap_betweenMetrics_thresh-%f.R.func.gii', parcelNum(iParcelNum), thresh)).cdata;

    % Concatenating the Hub Map
    hubMap = [hubMap_lh; hubMap_rh];

    % Searching for all the indices within the hub map that are nonzero (i.e.: a hub)
    idxNonZero = find(hubMap);

    % Getting the unique list of ROIs from the schaefer parcellation which are considered 'hubs' by our standards
    ROIlist = unique(parcelScheme(idxNonZero));
    
    % Creating the list
    list = {num2str(ROIlist', '%g,')};

    % Creating and writing the table for later use in python plotting
    T = table(parcelNum(iParcelNum), {list}, 'VariableNames', {'SchaeferParcellationNumber', 'ListOfHubROIs'});

    writetable(T, sprintf('/mfip/mfip1/arielle/PhDProject2/results/%s/hubLocations/hubROIsToPlot.csv', group), 'WriteMode', 'append')


end