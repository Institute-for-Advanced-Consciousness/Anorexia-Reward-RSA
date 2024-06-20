function [patterns_of_interest]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name, flags);

% Start Date: 4/24/24
% Contact: Nicco Reggente, Ph.D. (nicco@advancedconsciousness.org)

%This script will output a matrix with rows of features (i.e. voxel values)
%with columns of instances (i.e. trials / exemplars).
% This script is specific to a searchlight procedure to save time and
% memory.

%% Input Insights
%Inputs must take the following form and [var type]:

%subj [string] = Subject Number
%mask [vector] = Indices for active voxels (must be in same space as image)
%paths [struct] = Paths structure with relevant directory paths
%naming [struct] = Naming conventions for loading in subject data
%trials_of_interest [vector] = Indices for active trials to be loaded in
%from a subject's 4D file.

%% Create Searchlight
V=spm_vol([paths.masks '/' mask_name '.nii']);
vols=spm_read_vols(V);
[~ , ~, SL_Mask_Indices] = adj_sphere(vols);


%% Bring in Subject-Specific Data

%Load in Subject's 4D file
subj_4D_file=paths.data;
if ~exist(subj_4D_file)
    gunzip([subj_4D_file '.gz']); %Maybe was not gunzipped first.
end
V=spm_vol(subj_4D_file);
vols_of_interest=spm_read_vols(V);

%% Get Rid of motion outliers
active_trials=csvread(paths.active_trials);
bad_trials=find(active_trials==1);
if ~isempty(bad_trials)
vols_of_interest(:,:,:,bad_trials)=[];
end

clear vols V %Keep tabs on memory to run smoothly

%% SEARCHLIGHT PATTERN CREATION
disp('Starting To Loop Over Searchlights To Extract Data...')
for sl=1:size(SL_Mask_Indices,1)
    cur_sl_idx=find(SL_Mask_Indices(sl,:));
    mask=SL_Mask_Indices(sl,cur_sl_idx);
    for t=1:size(vols_of_interest,4)
        for m=1:length(mask)
            [a, b, c]=ind2sub(flags.default_size, mask(m));
            patterns_of_interest{sl}{t,m}=reshape_to_1(vols_of_interest(a,b,c,t));
        end
    end
    progress(sl,size(SL_Mask_Indices,1),10);
end
end