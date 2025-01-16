clc; clear all; close all; 

%% Loading paths
path = '/mfip/mfip1/arielle/PhDProject2/';
%parcellation = strsplit(input("What parcellation scheme and subject database do you want to do? (i.e.: schaefer 100 MICs HC) \n" , "s"));
    parcelName = 'schaefer'; %parcellation{1};
  %  parcelNum  = str2double(parcellation{2});
    groupName = 'MICs'; %parcellation{3};
  %  subjType = parcellation{4};

%% Loading Data - HC
% Loading Metric Info 
info_HC = readtable(fullfile(path, 'subject_lists', sprintf('subject_list_%s_HC_final.xlsx', groupName)));
sub_list_HC = info_HC.Subj_ID;
info_PX = readtable(fullfile(path, 'subject_lists', sprintf('subject_list_%s_PX_final.xlsx', groupName)));
sub_list_PX = info_PX.Subj_ID;

% Defining Parcellation Array
parcels = 100:100:900;
thresh  = [0.05:0.05:0.5];
subs    = {'HC', 'PX'}; 
%% Loading Data and Calculating Average Clustering Coefficient for HCs

for iSubType = 1:length(subs)

    if strcmp(subs{iSubType}, 'HC')
        sub_list = sub_list_HC; 
    elseif strcmp(subs{iSubType}, 'PX')
        sub_list = sub_list_PX;
    end

    for iParcel = 1:length(parcels)
    
        for ithresh = 1:length(thresh)
    
            for iSub = 1:length(sub_list)
                subData(:, iSub) = table2array(readtable(fullfile(path, 'data', groupName, 'clusteringCoefficient', sprintf('%s%d', parcelName, parcels(iParcel)), ...
                    sprintf('clusteringCoefficientResults_sub-%s_thresh-%f.csv', sub_list{iSub}, thresh(ithresh)))));
            end
            
            fprintf('Finished Loading %s for threshold %f at parcellation %d \n', subs{iSubType}, thresh(ithresh), parcels(iParcel))
            
            % transposing the data b/c i prefer doing math that way
            subData = subData';
    
            % Calculate the avg Clustering Coefficient for HC at the threshold thresh(ithresh) for parcellation schame parcels(iparcel) - aka computing the global clustering coefficient
            avgClusteringCoefficient = mean(mean(subData));
    
            T = table({sprintf('%s%d', parcelName, parcels(iParcel))}, thresh(ithresh), {subs{iSubType}}, avgClusteringCoefficient, ...
                'VariableNames', {'Parcellation', 'Threshold', 'Group', 'GlobalClusteringCoefficient'});
    
            writetable(T, fullfile(path, 'results', groupName, 'clusteringCoefficient', 'averageClusteringCoefficients.csv'), 'WriteMode', 'append')
    
            clear subData
        end
    
    end
end

