function Copy_of_BetweennessCentrality_BCT(parcelName, parcelNum, groupName)

addpath(genpath('/mfip/mfip1/arielle/software/BrainConnectivityToolbox'));
addpath(genpath('/mfip/mfip1/arielle/software/fdr_bh'));

%% For Visualization
addpath(genpath('/mfip/mfip1/arielle/software/BrainSpace'));
addpath(genpath('/mfip/mfip1/arielle/software/matlab_GIfTI'));
addpath(genpath('/mfip/mfip1/arielle/software/slanCM'));

[surf_lh, surf_rh] = load_conte69();
[mask_lh, mask_rh] = load_mask();

if nargin < 4
    parcellation = strsplit(input("What parcellation scheme and subject database do you want to do? (i.e.: schaefer 100 MICs HC) \n" , "s"));
    parcelName = parcellation{1};
    parcelNum  = str2double(parcellation{2});
    groupName = parcellation{3};
    subjType = parcellation{4};
end

% choose the number of threads for parfor
if ~isempty(getenv('NSLOTS'))
    thrds = maxNumCompThreads(str2num(getenv('NSLOTS'))) ;
    disp(thrds)
else
    thrds=10;
end

% open parfor if nescessary 
if isempty(gcp('nocreate'))
    parpool('Threads', thrds);
end

%% Getting user input for parcellation scheme to use
parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum)));
toVisualize = zeros(size(parcelScheme));

% Creating new directories
if ~exist(sprintf('/mfip/mfip1/arielle/PhDProject2/results/%s/%s%d/betweennessCentrality/', groupName, parcelName, parcelNum), 'dir')
    mkdir(sprintf('/mfip/mfip1/arielle/PhDProject2/results/%s/%s%d/betweennessCentrality/', groupName, parcelName, parcelNum))
end



%% Create Average Matrix
if strcmp(subjType, 'PX')
    info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_list_%s_final.xlsx', groupName));
    sub_list = info.Subj_ID;
elseif strcmp(subjType, 'HC')
    info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_list_%s_HC_final.xlsx', groupName));
    sub_list = info.Subj_ID;
end

for iSub = 1:size(info, 1)

  corrMatrix(:, :, iSub) = load(sprintf('/mfip/mfip1/arielle/PhDProject2/data/%s/correlationMatrixes/%s%d_correlationMatrixes/sub-%s_ses-01_surf-fsLR-32k_parc-%s%d_desc-corrMatrix.mat', ...
        groupName, parcelName, parcelNum, sub_list{iSub}, parcelName, parcelNum)).corrMatrix;

end

avgCorrMatrix = mean(corrMatrix,3);


%% Applying Graph Metrics

% setting thresholds we want to do
% thresh = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]; %, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1];
thresh = 0.1;
for ithresh = 1:length(thresh)

    % Thresholding the matrix
    threshMatrix =  threshold_proportional(avgCorrMatrix, thresh(ithresh));

    % computing the undirected connection-length matrix
    connectionLength = weight_conversion(threshMatrix, 'lengths');

    % compute the betweenness centrality : Nodes with high values of betweenness centrality participate in a large number of shortest paths. 
    % measures how often a node appears on the shortest paths between pairs of other nodes. 
    % calculated as: (the total number of shortest paths from node a to be) / (the number of the paths that pass through node c)
    % a way of defining hubness
    BetweennessCentrality=betweenness_wei(connectionLength)';

    % Saving results
    T = table(BetweennessCentrality', 'VariableNames', {'BetweennessCentrality'});
        writetable(T, sprintf('/mfip/mfip1/arielle/PhDProject2/results/%s/betweennessCentrality/%s%d_averageBetweennessCentralityResults_thresh-%f_group-%s.csv', ...
        groupName, parcelName, parcelNum, thresh(ithresh), subjType));
  
    % Visualization using brainspace - only visualizing nodes with significant p values
    toVisualize = zeros(size(parcelScheme));
    for idx = 1:size(avgCorrMatrix, 1)
         idxBig = find(parcelScheme == idx);
         toVisualize(idxBig) = BetweennessCentrality(idx);
    end
    
    % Normalizing Results
    toVisualize = toVisualize/max(toVisualize);

    % splitting the left and right hemispheres for plotting
    lh = toVisualize(1:size(toVisualize, 1)/2);
    rh = toVisualize((size(toVisualize, 1)/2)+1:end);
    
    % saving the gifti files
    savegifti(lh, rh, groupName, 'betweennessCentrality', parcelName, parcelNum, thresh, subjType)
end


end

%% Extra functions
function savegifti(mapToSaveLeft, mapToSaveRight, groupName, metric, parcelName, parcelNum, thresh, subjType)
    g = gifti(mapToSaveLeft);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/%s/%s_%s%d_thresh-%f_subjType-%s_averageMap.L.func.gii', groupName, metric, metric, parcelName, parcelNum, thresh, subjType))
    g = gifti(mapToSaveRight);
    save(g, sprintf('/mfip/mfip1/arielle/PhDProject2/maps/%s/%s/%s_%s%d_thresh-%f_subjType-%s_averageMap.R.func.gii', groupName, metric, metric, parcelName, parcelNum, thresh, subjType))
end