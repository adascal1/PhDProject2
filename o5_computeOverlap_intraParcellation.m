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
BC = table2array(readtable(fullfile(path, 'results', groupName, 'betweennessCentrality', sprintf('betweennessCentrality_%s%d_cohendResults_thresh-%f.csv', parcelName, parcelNum, thresh))));
EC = table2array(readtable(fullfile(path, 'results', groupName, 'eigenvectorCentrality', sprintf('eigenvectorCentrality_%s%d_cohendResults_thresh-%f.csv', parcelName, parcelNum, thresh))));

%% Finding where we have increases for all three metrics
% Finding just increases
DC_bin = DC > 0;
EB_bin = BC > 0;
EC_bin = EC > 0;

% putting all into one 3 x nParcel array 
bin = [DC_bin, EB_bin, EC_bin];
intersectAbove = sum(bin, 2);
consistentAbove = intersectAbove == 3;

%% Finding where we have decreases for all three metrics
% Finding just increases
DC_bin = DC < 0;
EB_bin = BC < 0;
EC_bin = EC < 0;

% putting all into one 3 x nParcel array 
bin = [DC_bin, EB_bin, EC_bin];
intersectBelow = sum(bin, 2);
consistentBelow = intersectBelow == 3;

toPercentage = [consistentAbove; consistentBelow];
consistencyPercentage = sum(toPercentage)/parcelNum;

consistencyAbovePercentage = sum(consistentAbove)/parcelNum;
consistencyBelowPercentage = sum(consistentBelow)/parcelNum;


T = table({sprintf('%s%d', parcelName, parcelNum)}, consistencyPercentage, consistencyAbovePercentage, consistencyBelowPercentage, 'VariableNames', {'Parcellation', 'TotalConsistencyPercentage', 'ConsistentlyAbovePercentage', 'ConsistentlyBelowPercentage'});
writetable(T, fullfile(path, 'results', 'MICs', 'consistency', 'consistencyOfMetrics_intraParcellation.csv'), 'WriteMode', 'append')
%% Creating Gifti Map

% creating gifti files - consistently above
toVisualize = zeros(size(parcelScheme));

for idx = 1:parcelNum

    idxBig = find(parcelScheme == idx); % finding the indexes in the parcel scheme associated with the parcel number 
    toVisualize(idxBig) = consistentAbove(idx); % placing the cohend at these index values

end

% splitting the left and right hemispheres for plotting
lh = toVisualize(1:size(toVisualize, 1)/2);
rh = toVisualize((size(toVisualize, 1)/2)+1:end);

% saving the gifti files
savegifti(lh, rh, groupName, parcelName, parcelNum, thresh, 'above');

% creating gifti files - consistently below
toVisualize = zeros(size(parcelScheme));

for idx = 1:parcelNum

    idxBig = find(parcelScheme == idx); % finding the indexes in the parcel scheme associated with the parcel number 
    toVisualize(idxBig) = consistentBelow(idx); % placing the cohend at these index values

end

% splitting the left and right hemispheres for plotting
lh = toVisualize(1:size(toVisualize, 1)/2);
rh = toVisualize((size(toVisualize, 1)/2)+1:end);

% saving the gifti files
savegifti(lh, rh, groupName, parcelName, parcelNum, thresh, 'below');





















%% Extra Functions
function savegifti(mapToSaveLeft, mapToSaveRight, groupName, parcelName, parcelNum, thresh, consistent)
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/%s%d_consistent-%s_thresh-%f.L.func.gii', groupName, parcelName, parcelNum, consistent, thresh))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/%s%d_consistent-%s_thresh-%f.R.func.gii', groupName, parcelName, parcelNum, consistent, thresh))
end
