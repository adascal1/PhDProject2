clc;
clear all;
close all;

addpath(genpath('/mfip/mfip1/arielle/software/matlab_GIfTI'));
addpath(genpath('/mfip/mfip1/arielle/software/BrainSpace'));

%% Getting user input for parcellation scheme to use
parcellation = strsplit(input("HC or PX, parcelName, SubjectGroup: \n" , "s"));
subjType = parcellation{1};
parcelName = parcellation{2};
subjGroup = parcellation{3}; %(MICS)

if strcmp(parcelName, 'schaefer')
    parcelNum  = 100:100:900;
end

%% Saving loading .gii file and saving as .mat file

% Setting paths / loading parcellations
micapath = '/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0/';
savepath = sprintf('/mfip/mfip1/arielle/PhDProject2/data/%s/correlationMatrixes', subjGroup);

if ~exist(savepath, 'dir')
    mkdir(savepath)
end

if strcmp(subjType, 'HC')
    info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_lists_%s_%s.xlsx', subjGroup, subjType));
    sub_list = info.Subj_ID;

elseif strcmp(subjType, 'PX')
      
    info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_lists_%s.xlsx', subjGroup));
    
    % removing bilateral temporal lobe cases and non defined cases
    info(strcmp(info.lateralization, "BL"), :) = [];
    info(strcmp(info.lateralization, "unc"), :) = [];
    sub_list = info.Subj_ID;

end

for iSub = 1:length(sub_list)

        if isfile(fullfile(micapath, sprintf('sub-%s', sub_list{iSub}), 'ses-01', 'func', 'desc-se_task-rest_acq-AP_bold', 'surf', sprintf('sub-%s_ses-01_surf-fsLR-32k_desc-timeseries_clean.shape.gii', sub_list{iSub})))
            subListNew{iSub} = sub_list{iSub};
        end
        
end

subListNew = subListNew(~cellfun('isempty', subListNew));

%keeping only the subjects with scans
[C, ia, ib] = intersect(info.Subj_ID, subListNew);
info = info(ia, :);
sub_list = info.Subj_ID;

writetable(info, sprintf('/mfip/mfip1/arielle/PhDProject2/subject_lists/subject_list_%s_%s_final.xlsx', subjGroup, subjType));

for iparcel = 1:length(parcelNum)

    parcelScheme = table2array(readtable(sprintf('/data/mica1/01_programs/micapipe-v0.2.0/parcellations/%s-%d_conte69.csv', parcelName, parcelNum(iparcel))));
    
    % indices for flipping purposes
    idx_right =  (1+parcelNum(iparcel)/2) : (parcelNum(iparcel));
    idx_left  =   1              : (parcelNum(iparcel)/2);
    
    
    %% Saving corr matrices and performing flipping
    for iSub = 1:size(info, 1)
    
         if isfile(fullfile(savepath, sprintf('%s%d_correlationMatrixes', parcelName, parcelNum(iparcel)), sprintf('sub-%s_ses-01_surf-fsLR-32k_parc-%s%d_desc-corrMatrix.mat', sub_list{iSub}, parcelName, parcelNum(iparcel))))
             fprintf('Subject %s already computed, next! \n', sub_list{iSub})
         else

            % Load the timeseries
            ts = fullfile(micapath, sprintf('sub-%s', sub_list{iSub}), 'ses-01', 'func', 'desc-se_task-rest_acq-AP_bold', 'surf', sprintf('sub-%s_ses-01_surf-fsLR-32k_desc-timeseries_clean.shape.gii', sub_list{iSub}));
            ts = gifti(ts).cdata;
        
            % parcellate the time series
            ts_parcel = full2parcel(ts, parcelScheme)';
        
            % flip rTLEs to all seizure foci are on the left hemisphere
            if strcmp(subjType, 'PX')
                if strcmp(info{iSub, 'lateralization'}, 'R')
                    ts_parcel = ts_parcel(:, [idx_right, idx_left]);
                end
            end
        
            % creating the correlation matrix
            corrMatrix = corr(ts_parcel);
        
            % saving the correlation matrix as a .mat file
            filename = fullfile(savepath, sprintf('%s%d_correlationMatrixes', parcelName, parcelNum(iparcel)), sprintf('sub-%s_ses-01_surf-fsLR-32k_parc-%s%d_desc-corrMatrix.mat', sub_list{iSub}, parcelName, parcelNum(iparcel)));
            save(filename, "corrMatrix");
        
            fprintf('finished subject %s \n', sub_list{iSub})

        end
    end
    
    fprintf('finished parcel %d \n', parcelNum(iparcel))
end
