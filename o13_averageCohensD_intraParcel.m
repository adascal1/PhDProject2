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

parcelNums = {'100', '300', '600', '900'};

%% Loading the Maps

for iparcelNum = 1:length(parcelNums)

    lh_DC = gifti(fullfile(path, 'maps', groupName, 'degree', sprintf('degree_%s%s_thresh-%f_cohend.L.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;
    rh_DC = gifti(fullfile(path, 'maps', groupName, 'degree', sprintf('degree_%s%s_thresh-%f_cohend.R.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;

    lh_BC = gifti(fullfile(path, 'maps', groupName, 'betweennessCentrality', sprintf('betweennessCentrality_%s%s_thresh-%f_cohend.L.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;
    rh_BC = gifti(fullfile(path, 'maps', groupName, 'betweennessCentrality', sprintf('betweennessCentrality_%s%s_thresh-%f_cohend.R.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;

    lh_EC = gifti(fullfile(path, 'maps', groupName, 'eigenvectorCentrality', sprintf('eigenvectorCentrality_%s%s_thresh-%f_cohend.L.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;
    rh_EC = gifti(fullfile(path, 'maps', groupName, 'eigenvectorCentrality', sprintf('eigenvectorCentrality_%s%s_thresh-%f_cohend.R.func.gii', parcelName, parcelNums{iparcelNum}, thresh))).cdata;

    lh = [lh_DC, lh_BC, lh_EC];
    rh = [rh_DC, rh_BC, rh_EC];
    
    mean_lh = mean(lh, 2);
    mean_rh = mean(rh, 2);

    savegifti(mean_lh, mean_rh, groupName, parcelName, parcelNums{iparcelNum}, thresh);

    
end




%% Extra Functions
function savegifti(mapToSaveLeft, mapToSaveRight, groupName, parcelName, parcelNum, thresh)
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/%s%s_avgCohenD_intraParcel_thresh-%f.L.func.gii', groupName, parcelName, parcelNum, thresh))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/%s%s_avgCohenD_intraParcel_thresh-%f.R.func.gii', groupName, parcelName, parcelNum, thresh))
end