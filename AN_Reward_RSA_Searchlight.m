function []=AN_Reward_RSA_Searchlight(subj,mask_name,analysis_type)
% Start Date: 4/24/24
% Contact: Nicco Reggente, Ph.D. (nicco@advancedconsciousness.org)

%Workspace Running
%AN_Reward_RSA_Searchlight('4011','Amygdala_Basolateral',1)

%% Input Insights
%Inputs must take the following form and [var type]:

%data_types [cell array of numbers] = Which data types (max 2) will be used
%to run the RSA comparison of within vs. between. (e.g. {1,2}

%analysis_type [numeric]
% 1 = Reward Responses during Neutral States (All Subjs --> AN vs. Control)
%     "Reward_Only_Trials" & "LSSout_MNI_NEU"
%
% 2 = Reward Responses during Anxious States (All Subjs --> AN vs. Control)
%     "Reward_Only_Trials" & "LSSout_MNI_ANX"
%
% 3 = Anxiety vs. Neutral (All Subjs --> AN Clinical Outcomes)
%     "Reward_Only_Trials" & "LSSout_MNI_ANX" vs. "LSSout_MNI_NEU"
%
% 4 = Anxiety Word Presentation (All Subjs --> AN vs. Control)
%     "Word_After_Reward" & "LSSout_MNI_ANX"
%
% 5 = Non-Reward Responses during Neutral States (All Subjs --> AN vs. Control)
%     "Non_Reward_Only_Trials" & "LSSout_MNI_NEU"
%
% 6 = Non-Reward Responses during Anxious States (All Subjs --> AN vs. Control)
%     "Non_Reward_Only_Trials" & "LSSout_MNI_ANX"
%
% 7 = Anxiety vs. Neutral (All Subjs --> AN Clinical Outcomes)
%     "Non_Reward_Only_Trials" & "LSSout_MNI_ANX" vs. "LSSout_MNI_NEU"
%
% 8 = Neutral Reward vs. Neutral Non Reward (All Subjs --> AN vs. Control)
%     "Non_Reward_Only_Trials & LSSout_MNI_NEU" vs. "Reward_Only_Trials & LSSout_MNI_NEU"


disp(['Running Subject....' subj])
%% Hard Code
flags.default_size=[91 109 91];

%% Leverage Naming Conventions
naming.pre_subj='sub_'; %The standard prefix on files. Include underscores if they are used.
naming.subj_folder='sub-';

%% Define Paths
paths.top='D:\PROJECTS\AN_Reward_RSA\';
paths.reference=[paths.top 'Reference'];
paths.masks=[paths.top 'Masks/'];
paths.ml='D:\PROJECTS\MATLAB_PATH'; addpath(genpath(paths.ml)); %add toolboxes to the path
paths.save=[paths.top 'Subject_Results/Searchlights/' mask_name '/'];
if ~exist(paths.save)
    mkdir(paths.save)
end

%% Initiate Analysis
switch analysis_type
    case 1

        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_reward_motion/sub-' num2str(subj) '.csv'];

        [patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);

        %% Create Within Maps
        disp('Creating Within Reward (Neutral) Patterns...')
        for sl=1:size(patterns,2)% Loop over searchlights
            A=tril(corr(cell2mat(patterns{sl})'),-1);%Within each searchlight sphere, correlate all the exmemplars
            temp_r{sl}=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
            progress(sl,size(patterns,2),10);
        end
        Within_Reward_Neutral_Values=cell2mat(temp_r);
        %     delete(ppm);
        clear temp_r

        %% Write Outputs

        V=spm_vol([paths.masks '/' mask_name '.nii']);
        temp_vols=spm_read_vols(V);
        active_mask_idx=find(temp_vols);
        V.dt=[16 0];

        %Write Within Me Maps
        vols=temp_vols;
        vols(active_mask_idx)=0;%First clear it out.
        vols(active_mask_idx)=Within_Reward_Neutral_Values; %assign R values to the active voxels.
        V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Reward_Neutral.nii'];
        spm_write_vol(V,vols);

    case 2

        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_ANX.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/anxiety_reward_motion/sub-' num2str(subj) '.csv'];

        [patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);

        %% Create Within Maps
        disp('Creating Within Reward (Anxiety) Patterns...')
        for sl=1:size(patterns,2)% Loop over searchlights
            A=tril(corr(cell2mat(patterns{sl})'),-1);%Within each searchlight sphere, correlate all the exmemplars
            temp_r{sl}=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
            progress(sl,size(patterns,2),10);
        end
        Within_Reward_Anxiety_Values=cell2mat(temp_r);
        %     delete(ppm);
        clear temp_r

        %% Write Outputs
        V=spm_vol([paths.masks '/' mask_name '.nii']);
        temp_vols=spm_read_vols(V);
        active_mask_idx=find(temp_vols);
        V.dt=[16 0];

        %Write Within Me Maps
        vols=temp_vols;
        vols(active_mask_idx)=0;%First clear it out.
        vols(active_mask_idx)=Within_Reward_Anxiety_Values; %assign R values to the active voxels.
        V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Reward_Anxiety.nii'];
        spm_write_vol(V,vols);


    case 3
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_ANX.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/anxiety_reward_motion/sub-' num2str(subj) '.csv'];
        [anxiety_patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);
        paths.data=[paths.top 'Data/Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_reward_motion/sub-' num2str(subj) '.csv'];
        [neutral_patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);
        %% Create Across Maps
        disp('Creating Across Reward (Anxiety vs. Neutral) Patterns...')
        anxiety_patterns=anxiety_patterns';
        neutral_patterns=neutral_patterns';

        for sl=1:size(anxiety_patterns,1)% Loop over searchlights
            A=size(anxiety_patterns{sl},1);
            B=size(neutral_patterns{sl},1);
            C=zeros(1,A*B);
            for a=1:A
                for b=1:B
                    C(1,(a*b))=corr(cell2mat(anxiety_patterns{sl}(a,:)'),cell2mat(neutral_patterns{sl}(b,:)'));
                end
            end
            temp_r{sl}=mean(atanh(C(find(C))));%Grab the non-zero values and do an r-to-z transform.
            progress(sl,size(anxiety_patterns,1),10);
        end

        Across_Reward_Values=cell2mat(temp_r);
        clear temp_r

        %% Write Outputs
        V=spm_vol([paths.masks '/' mask_name '.nii']);
        temp_vols=spm_read_vols(V);
        active_mask_idx=find(temp_vols);
        V.dt=[16 0];

        %Write Within Me Maps
        vols=temp_vols;
        vols(active_mask_idx)=0;%First clear it out.
        vols(active_mask_idx)=Across_Reward_Values; %assign R values to the active voxels.
        V.fname=[paths.save subj '_' mask_name '_Searchlight_Across_Reward_Anxiety_vs_Neutral.nii'];
        spm_write_vol(V,vols);


    case 4

        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Word_After_Reward/' naming.subj_folder num2str(subj) '/LSSout_MNI_ANX.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/anx_after_reward/sub-' num2str(subj) '_REWARD_ANX.csv'];

        [patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);

        %% Create Within Maps
        disp('Creating Within Anxiety Word Patterns...')
        for sl=1:size(patterns,2)% Loop over searchlights
            A=tril(corr(cell2mat(patterns{sl})'),-1);%Within each searchlight sphere, correlate all the exmemplars
            temp_r{sl}=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
            progress(sl,size(patterns,2),10);
        end
        Within_Anxiety_Word_Values=cell2mat(temp_r);
        %     delete(ppm);
        clear temp_r

        %% Write Outputs
        V=spm_vol([paths.masks '/' mask_name '.nii']);
        temp_vols=spm_read_vols(V);
        active_mask_idx=find(temp_vols);
        V.dt=[16 0];

        %Write Within Me Maps
        vols=temp_vols;
        vols(active_mask_idx)=0;%First clear it out.
        vols(active_mask_idx)=Within_Anxiety_Word_Values; %assign R values to the active voxels.
        V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Anxiety_Word.nii'];
        spm_write_vol(V,vols);

    case 5
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Non_Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_non-reward_motion/sub-' num2str(subj) '.csv'];

        [patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);

        %% Create Within Maps
        disp('Creating Within Non Reward (Neutral) Patterns...')
        for sl=1:size(patterns,2)% Loop over searchlights
            A=tril(corr(cell2mat(patterns{sl})'),-1);%Within each searchlight sphere, correlate all the exmemplars
            temp_r{sl}=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
            progress(sl,size(patterns,2),10);
        end
        Within_NonReward_Neutral_Values=cell2mat(temp_r);
        %     delete(ppm);
        clear temp_r

        %% Write Outputs

        V=spm_vol([paths.masks '/' mask_name '.nii']);
        temp_vols=spm_read_vols(V);
        active_mask_idx=find(temp_vols);
        V.dt=[16 0];

        %Write Within Me Maps
        vols=temp_vols;
        vols(active_mask_idx)=0;%First clear it out.
        vols(active_mask_idx)=Within_NonReward_Neutral_Values; %assign R values to the active voxels.
        V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_NonReward_Neutral.nii'];
        spm_write_vol(V,vols);
    case 6
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Non_Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_ANX.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/anxiety_non-reward_motion/sub-' num2str(subj) '.csv'];

        [patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);

        %% Create Within Maps
        disp('Creating Within Non Reward (Anxiety) Patterns...')
        for sl=1:size(patterns,2)% Loop over searchlights
            A=tril(corr(cell2mat(patterns{sl})'),-1);%Within each searchlight sphere, correlate all the exmemplars
            temp_r{sl}=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
            progress(sl,size(patterns,2),10);
        end
        Within_NonReward_Anxiety_Values=cell2mat(temp_r);
        %     delete(ppm);
        clear temp_r

        %% Write Outputs
        V=spm_vol([paths.masks '/' mask_name '.nii']);
        temp_vols=spm_read_vols(V);
        active_mask_idx=find(temp_vols);
        V.dt=[16 0];

        %Write Within Me Maps
        vols=temp_vols;
        vols(active_mask_idx)=0;%First clear it out.
        vols(active_mask_idx)=Within_NonReward_Anxiety_Values; %assign R values to the active voxels.
        V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_NonReward_Anxiety.nii'];
        spm_write_vol(V,vols);

    case 7
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Non_Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_ANX.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/anxiety_non-reward_motion/sub-' num2str(subj) '.csv'];
        [anxiety_patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);
        paths.data=[paths.top 'Data/Non_Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_non-reward_motion/sub-' num2str(subj) '.csv'];

        [neutral_patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);
        %% Create Across Maps
        disp('Creating Across Non Reward (Anxiety vs. Neutral) Patterns...')
        anxiety_patterns=anxiety_patterns';
        neutral_patterns=neutral_patterns';

        for sl=1:size(anxiety_patterns,1)% Loop over searchlights
            A=size(anxiety_patterns{sl},1);
            B=size(neutral_patterns{sl},1);
            C=zeros(1,A*B);
            for a=1:A
                for b=1:B
                    C(1,(a*b))=corr(cell2mat(anxiety_patterns{sl}(a,:)'),cell2mat(neutral_patterns{sl}(b,:)'));
                end
            end
            temp_r{sl}=mean(atanh(C(find(C))));%Grab the non-zero values and do an r-to-z transform.
            progress(sl,size(anxiety_patterns,1),10);
        end

        Across_NonReward_Values=cell2mat(temp_r);
        clear temp_r

        %% Write Outputs
        V=spm_vol([paths.masks '/' mask_name '.nii']);
        temp_vols=spm_read_vols(V);
        active_mask_idx=find(temp_vols);
        V.dt=[16 0];

        %Write Within Me Maps
        vols=temp_vols;
        vols(active_mask_idx)=0;%First clear it out.
        vols(active_mask_idx)=Across_NonReward_Values; %assign R values to the active voxels.
        V.fname=[paths.save subj '_' mask_name '_Searchlight_Across_NonReward_Anxiety_vs_Neutral.nii'];
        spm_write_vol(V,vols);

    case 8
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Non_Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_non-reward_motion/sub-' num2str(subj) '.csv'];

        [non_reward_neutral_patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);
        paths.data=[paths.top 'Data/Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_reward_motion/sub-' num2str(subj) '.csv'];

        [reward_neutral_patterns]=AN_Reward_RSA_create_SL_patterns(subj,paths,mask_name,flags);
        %% Create Across Maps
        disp('Creating Across Non Reward (Anxiety vs. Neutral) Patterns...')
        non_reward_neutral_patterns=non_reward_neutral_patterns';
        reward_neutral_patterns=reward_neutral_patterns';

        for sl=1:size(non_reward_neutral_patterns,1)% Loop over searchlights
            A=size(non_reward_neutral_patterns{sl},1);
            B=size(reward_neutral_patterns{sl},1);
            C=zeros(1,A*B);
            for a=1:A
                for b=1:B
                    C(1,(a*b))=corr(cell2mat(non_reward_neutral_patterns{sl}(a,:)'),cell2mat(reward_neutral_patterns{sl}(b,:)'));
                end
            end
            temp_r{sl}=mean(atanh(C(find(C))));%Grab the non-zero values and do an r-to-z transform.
            progress(sl,size(non_reward_neutral_patterns,1),10);
        end

        Across_Neutral_Values=cell2mat(temp_r);
        clear temp_r

        %% Write Outputs
        V=spm_vol([paths.masks '/' mask_name '.nii']);
        temp_vols=spm_read_vols(V);
        active_mask_idx=find(temp_vols);
        V.dt=[16 0];

        %Write Within Me Maps
        vols=temp_vols;
        vols(active_mask_idx)=0;%First clear it out.
        vols(active_mask_idx)=Across_Neutral_Values; %assign R values to the active voxels.
        V.fname=[paths.save subj '_' mask_name '_Searchlight_Across_NonReward_Neutral_vs_Reward_Neutral.nii'];
        spm_write_vol(V,vols);
end
end

%
%
% %% Create Across Maps
% disp('Creating Across Patterns...')
% for p=1:(size(patterns,2)-1);%Loop over each type of pair.
%     Temp1=patterns{p};
%     Temp2=patterns{p+1};
%     %     ppm=ParforProgressbar(size(Temp1,2),'showWorkerProgress',true,'progressBarUpdatePeriod',4,'title','Across Pattern Building');
%     for sl=1:size(patterns{p},2)% Loop over searchlights
%         A=size(Temp1{sl},1);
%         B=size(Temp2{sl},1);
%         C=zeros(1,A*B);
%         for a=1:A
%             for b=1:B
%                 C(1,(a*b))=corr(cell2mat(Temp1{sl}(a,:))',cell2mat(Temp2{sl}(b,:))');
%             end
%         end
%         temp_r{sl}=mean(atanh(C));%Grab the non-zero values and do an r-to-z transform.
%         %         ppm.increment()
%         progress(sl,size(patterns{p},2),20);
%     end
%     %     delete(ppm);
%     Across_R_Values{p}=cell2mat(temp_r);
%     clear temp_r
% end
%
%
% %% Create Within Maps
% disp('Creating Within Patterns...')
%
% for p=1:size(patterns,2)%Loop over each type of pattern.
%     %     ppm=ParforProgressbar(size(patterns{p},2),'showWorkerProgress',true,'progressBarUpdatePeriod',4,'title','Within Pattern Building');
%     for sl=1:size(patterns{p},2)% Loop over searchlights
%         A=tril(corr(cell2mat(patterns{p}{sl})'),-1);%Within each searchlight sphere, correlate all the exmemplars
%         temp_r{sl}=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
%         %         ppm.increment();
%         progress(sl,size(patterns{p},2),20);
%     end
%     Within_R_Values{p}=cell2mat(temp_r);
%     %     delete(ppm);
%     clear temp_r
%
% end
%
%
% %% Create Within vs. Across Comparisons
% %Statistically needs to be done at the group level, but can create the
% %subtraction maps here.
% across_map=Across_R_Values{1};
% within_anxiety_map=Within_R_Values{1};
% within_neutral_map=Within_R_Values{2};
%
% within_anxiety_map_minus_within_neutral_map=within_anxiety_map-within_neutral_map;
% within_neutral_map_minus_within_anxiety_map=within_neutral_map-within_anxiety_map;
%
% across_minus_within_anxiety_map=across_map-within_anxiety_map;
% across_minus_within_neutral_map=across_map-within_neutral_map;
%
% within_anxiety_minus_across=within_anxiety_map-across_map;
% within_neutral_minus_across=within_neutral_map-across_map;
%
%
% %% Write Outputs
%
% V=spm_vol([paths.masks '/' mask_name '.nii']);
% temp_vols=spm_read_vols(V);
% active_mask_idx=find(temp_vols);
% V.dt=[16 0];
%
% %Write Within Me Maps
% vols=temp_vols;
% vols(active_mask_idx)=within_anxiety_map; %assign R values to the active voxels.
% V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Anxiety.nii'];
% spm_write_vol(V,vols);
%
% %Write Within NotMe Maps
% vols=temp_vols;
% vols(active_mask_idx)=within_neutral_map; %assign R values to the active voxels.
% V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Neutral.nii'];
% spm_write_vol(V,vols);
%
% %Write Across Maps
% vols=temp_vols;
% vols(active_mask_idx)=across_map; %assign R values to the active voxels.
% V.fname=[paths.save subj '_' mask_name '_Searchlight_Across.nii'];
% spm_write_vol(V,vols);
%
% %Write Subtraction Maps
% vols=temp_vols;
% vols(active_mask_idx)=within_anxiety_map_minus_within_neutral_map; %assign R values to the active voxels.
% V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Anxiety_Minus_Within_Neutral.nii'];
% spm_write_vol(V,vols);
%
% vols=temp_vols;
% vols(active_mask_idx)=within_neutral_map_minus_within_anxiety_map; %assign R values to the active voxels.
% V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Neutral_Minus_Within_Anxiety.nii'];
% spm_write_vol(V,vols);
%
% vols=temp_vols;
% vols(active_mask_idx)=across_minus_within_anxiety_map; %assign R values to the active voxels.
% V.fname=[paths.save subj '_' mask_name '_Searchlight_Across_Minus_Within_Anxiety.nii'];
% spm_write_vol(V,vols);
%
% vols=temp_vols;
% vols(active_mask_idx)=across_minus_within_neutral_map; %assign R values to the active voxels.
% V.fname=[paths.save subj '_' mask_name '_Searchlight_Across_Minus_Within_Neutral.nii'];
% spm_write_vol(V,vols);
%
% vols=temp_vols;
% vols(active_mask_idx)=within_anxiety_minus_across; %assign R values to the active voxels.
% V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Anxiety_Minus_Across.nii'];
% spm_write_vol(V,vols);
%
% vols=temp_vols;
% vols(active_mask_idx)=within_neutral_minus_across; %assign R values to the active voxels.
% V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Neutral_Minus_Across.nii'];
% spm_write_vol(V,vols);
%
% disp(['Done with Subject ' subj '. That Took: ' num2str(toc/60) ' minutes.'])
%
% end