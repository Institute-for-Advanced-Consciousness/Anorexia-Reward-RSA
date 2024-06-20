function [patterns_of_interest]=create_ROI_patterns(paths, mask_name, flags)
% Start Date: 4/23/24
% Contact: Nicco Reggente, Ph.D. (nicco@advancedconsciousness.org)

%This script will output a matrix with rows of features (i.e. voxel values)
%with columns of instances (i.e. trials / exemplars).
% This script is specific to an ROI procedure

%% Input Insights
%Inputs must take the following form and [var type]:

%subj [string] = Subject Number
%mask [vector] = Indices for active voxels (must be in same space as image)
%paths [struct] = Paths structure with relevant directory paths
%naming [struct] = Naming conventions for loading in subject data
%trials_of_interest [vector] = Indices for active trials to be loaded in
%from a subject's 4D file.

%% Bring in Subject-Specific Data

%Load in Subject's 4D file
subj_4D_file=paths.data;
if ~exist(subj_4D_file)
    gunzip([subj_4D_file '.gz']); %Maybe was not gunzipped first.
end
V=spm_vol(subj_4D_file);
vols=spm_read_vols(V);

%% Mask the Data with voxels in a mask
V=spm_vol([paths.masks mask_name '.nii']);
mask_vols=spm_read_vols(V);
mask_indices=find(mask_vols);

%% Get Rid of motion outliers
active_trials=csvread(paths.active_trials);
bad_trials=find(active_trials==1);
if ~isempty(bad_trials)
vols(:,:,:,bad_trials)=[];
end

disp('Creating Patterns for ROI...')
for t=1:size(vols,4)
    for m=1:length(mask_indices)
        [a, b, c]=ind2sub(flags.default_size, mask_indices(m));
        patterns_of_interest{m,t}=reshape_to_1(vols(a,b,c,t));
    end
    progress(t,size(vols,4),10);
end
patterns_of_interest=cell2mat(patterns_of_interest);
end