clc; clear all; close all;

%% Loading Data
path = '/mfip/mfip1/arielle/PhDProject2/';
addpath(genpath('/mfip/mfip1/arielle/software/matlab_GIfTI'));
addpath(genpath('/mfip/mfip1/arielle/software/fdr_bh'));

% User input for analysis
parcellation = strsplit(input("What parcellation scheme and subject database do you want to do and metric? (i.e.: schaefer 100 MICs betweennessCentrality) \n" , "s"));
    parcelName  = parcellation{1};
    parcelNum   = str2double(parcellation{2});
    groupName   = parcellation{3};
    metric      = parcellation{4};
    thresh      = 0.1 ;%str2double(parcellation{5});

% Loading Metric Info 
info_HC = readtable(fullfile(path, 'subject_lists', sprintf('subject_list_%s_HC_final.xlsx', groupName)));
sub_list_HC = info_HC.Subj_ID;
info_PX = readtable(fullfile(path, 'subject_lists', sprintf('subject_list_%s_PX_final.xlsx', groupName)));
sub_list_PX = info_PX.Subj_ID;

for iSub = 1:length(sub_list_HC)
   HC(:, iSub) = table2array(readtable(fullfile(path, 'data', groupName, metric, sprintf('%s%d', parcelName, parcelNum), ...
       sprintf('%sResults_sub-%s_thresh-%f.csv', metric, sub_list_HC{iSub}, thresh))));
end

for iSub = 1:length(sub_list_PX)
   PX(:, iSub) = table2array(readtable(fullfile(path, 'data', groupName, metric, sprintf('%s%d', parcelName, parcelNum), ...
       sprintf('%sResults_sub-%s_thresh-%f.csv', metric, sub_list_PX{iSub}, thresh))));
end

% Transposing Data - creating an array of size (# of subjects in the rows, and # of parcels in columns)
HC = HC'; 
PX = PX'; 

% Loading parcellation
parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum)));


%% Computing Effect Size
% Creating new directories
if ~exist(sprintf('/mfip/mfip1/arielle/PhDProject2/results/%s/%s/', groupName, metric), 'dir')
    mkdir(sprintf('/mfip/mfip1/arielle/PhDProject2/results/%s/%s/', groupName, metric))
end

if ~exist(sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/%s/', groupName, metric), 'dir')
    mkdir(sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/%s/', groupName, metric))
end

% Effect size computation : d = mean(patients) - mean(healthy) / pooled standard deviation

% Calculating Mean metric at each parcel
mean_HC = mean(HC);
mean_PX = mean(PX);

% Calculating Pooled Standard Deviation at each parcel
pooledstd = sqrt(  ( (size(sub_list_PX, 1)-1).*std(PX).^2 + (size(sub_list_HC, 1)-1).*std(HC).^2 )  ./  ( (length(sub_list_HC) + length(sub_list_PX))-2 )  ); % calculating the pooled standard deviation for each vertex

% positive cohend = PX has a higher mean than HC 
% negative cohend = PX had a lower mean than HC
cohend = (mean_PX - mean_HC) ./ pooledstd;

%% Compute t-statistic between HC and PX : equal or unequal sample sizes with similar variance

% tstatistic
tstatistic = (mean_HC - mean_PX) ./ (pooledstd .* sqrt(1/length(sub_list_HC) + 1/length(sub_list_PX))); 

% degree of freedom
df = (length(sub_list_HC) + length(sub_list_PX))-2;

% pvalues
pValue = 2 * (1 - tcdf(abs(tstatistic), df));

% fdr Correction
[h, ~, ~, adj_p] = fdr_bh(pValue);


%% Saving table with cohend and maps on the surface
% saving data
fname = fullfile(path, 'results', groupName, metric, sprintf('%s_%s%d_cohendResults_thresh-%f.csv', metric, parcelName, parcelNum, thresh));
writetable(table(cohend', tstatistic', pValue', adj_p', h', 'VariableNames', {'cohenD', 'tstatistic', 'pVal', 'fdrCorrected', 'significance'}), fname)

% creating gifti files
toVisualize = zeros(size(parcelScheme));

for idx = 1:parcelNum

    idxBig = find(parcelScheme == idx); % finding the indexes in the parcel scheme associated with the parcel number 
    toVisualize(idxBig) = cohend(idx); % placing the cohend at these index values

end

% splitting the left and right hemispheres for plotting
lh = toVisualize(1:size(toVisualize, 1)/2);
rh = toVisualize((size(toVisualize, 1)/2)+1:end);

% saving the gifti files
savegifti(lh, rh, groupName, metric, parcelName, parcelNum, thresh)

%% Saving cohend map but only significant pvals
% creating gifti files
toVisualize = zeros(size(parcelScheme));

for idx = 1:parcelNum
    
    if adj_p(idx) < 0.05
        idxBig = find(parcelScheme == idx); % finding the indexes in the parcel scheme associated with the parcel number 
        toVisualize(idxBig) = cohend(idx); % placing the cohend at these index values
    end

end

% splitting the left and right hemispheres for plotting
lh = toVisualize(1:size(toVisualize, 1)/2);
rh = toVisualize((size(toVisualize, 1)/2)+1:end);

% saving the gifti files
savegifti_sig(lh, rh, groupName, metric, parcelName, parcelNum, thresh)

%% Extra Functions
function savegifti(mapToSaveLeft, mapToSaveRight, groupName, metric, parcelName, parcelNum, thresh)
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/%s/%s_%s%d_thresh-%f_cohend.L.func.gii', groupName, metric, metric, parcelName, parcelNum, thresh))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/%s/%s_%s%d_thresh-%f_cohend.R.func.gii', groupName, metric, metric, parcelName, parcelNum, thresh))
end

function savegifti_sig(mapToSaveLeft, mapToSaveRight, groupName, metric, parcelName, parcelNum, thresh)
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/%s/%s_%s%d_thresh-%f_cohend_sigROIs.L.func.gii', groupName, metric, metric, parcelName, parcelNum, thresh))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/%s/%s_%s%d_thresh-%f_cohend_sigROIs.R.func.gii', groupName, metric, metric, parcelName, parcelNum, thresh))
end
