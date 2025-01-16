clc; clear all; close all; 

path = '/mfip/mfip1/arielle/PhDProject2/';
addpath(genpath('/mfip/mfip1/arielle/software/matlab_GIfTI'));

% User input for analysis
%parcellation = strsplit(input("What parcellation scheme and subject database do you want to do and metric? (i.e.: schaefer 100 MICs 0.1) \n" , "s"));
    parcelName  = 'schaefer'; %parcellation{1};
    % parcelNum   = 100; %str2double(parcellation{2});
    groupName   = 'MICs'; %parcellation{3};
    thresh      = 0.1; %str2double(parcellation{4});

    parcellations = [100, 300, 600, 900];
    metrics       = {'betweennessCentrality', 'degree', 'eigenvectorCentrality', 'participationCoefficient'};

k = 1;
% loading all the parcellations
for iparcel = 1:length(parcellations)

    % loading all the metrics
    for imetric = 1:length(metrics)
        
        map_lh = gifti(fullfile(path, 'maps', groupName, metrics{imetric}, sprintf('%s_%s%d_thresh-%f_cohend.L.func.gii', metrics{imetric}, parcelName, parcellations(iparcel), thresh))).cdata;
        map_rh = gifti(fullfile(path, 'maps', groupName, metrics{imetric}, sprintf('%s_%s%d_thresh-%f_cohend.R.func.gii', metrics{imetric}, parcelName, parcellations(iparcel), thresh))).cdata;
        tmp = [map_lh; map_rh];

        mapArray(:, k) = tmp;
        k = k+1;

        load_order(k) = {sprintf('%s%d-%s', parcelName, parcellations(iparcel), metrics{imetric})};

    end

end
   

[correlationMatrix_byParcel, pvals] = corr(mapArray);



correlationMatrix_byMetric = corr(mapArray);

g = gifti(correlationMatrix_byMetric);
save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/consistency/correlationMatrixBetweenMaps_thresh-%f.L.func.gii', groupName, thresh))

   