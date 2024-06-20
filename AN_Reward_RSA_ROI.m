function []=AN_Reward_RSA_ROI(subj,mask_name,analysis_type)
% Start Date: 4/23/24
% Contact: Nicco Reggente, Ph.D. (nicco@advancedconsciousness.org)

%Workspace Running
%AN_RSA_ROI('4011','Amygdala_Basolateral',1)

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
paths.save=[paths.top 'Subject_Results/ROIs/'];

%% Initiate Analysis
switch analysis_type
    case 1 %Neutral Reward
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_reward_motion/sub-' num2str(subj) '.csv'];
        [patterns]=create_ROI_patterns(paths,mask_name,flags);

        %% Create Within Values
        disp('Creating Within Patterns For Neutral Reward Trials...')
        A=tril(corr(patterns),-1);%Within ROI, correlate all the exmemplars
        temp_r=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
        Within_R_Values=temp_r;

        %% Write Outputs
        header={'Within_Neutral_Reward_Trials'};
        data=Within_R_Values;
        save_file=[paths.save subj '_Within_Neutral_Reward_' mask_name '.txt'];
        save_data_with_headers(header,data,save_file);

    case 2 %Anxiety Reward
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_ANX.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/anxiety_reward_motion/sub-' num2str(subj) '.csv'];
        [patterns]=create_ROI_patterns(paths,mask_name,flags);

        %% Create Within Values
        disp('Creating Within Patterns For Anxiety Reward Trials...')
        A=tril(corr(patterns),-1);%Within ROI, correlate all the exmemplars
        temp_r=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
        Within_R_Values=temp_r;

        %% Write Outputs
        header={'Within_Anxiety_Reward_Trials'};
        data=Within_R_Values;
        save_file=[paths.save subj '_Within_Anxiety_Reward_' mask_name '.txt'];
        save_data_with_headers(header,data,save_file);

    case 3 %Neutral & Anxiety Reward
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_ANX.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/anxiety_reward_motion/sub-' num2str(subj) '.csv'];
        [anxiety_patterns]=create_ROI_patterns(paths,mask_name,flags);
        paths.data=[paths.top 'Data/Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_reward_motion/sub-' num2str(subj) '.csv'];
        [neutral_patterns]=create_ROI_patterns(paths,mask_name,flags);

        %% Create Across Values
        disp('Creating Across Patterns...')

        anxiety_patterns=anxiety_patterns';
        neutral_patterns=neutral_patterns';
        A=size(anxiety_patterns,1);
        B=size(neutral_patterns,1);
        C=zeros(1,A*B);
        for a=1:A
            for b=1:B
                C(1,(a*b))=corr(anxiety_patterns(a,:)',neutral_patterns(b,:)');
            end
            progress(a,A,20);
        end
        temp_r=mean(atanh(C(find(C))));%Grab the non-zero values and do an r-to-z transform.
        Across_R_Values=temp_r;
        clear temp_r

        %% Write Outputs
        header={'Across_Reward_Anxiety_and_Neutral_Trials'};
        data=Across_R_Values;
        save_file=[paths.save subj '_Across_Reward_Anxiety_and_Neutral_' mask_name '.txt'];
        save_data_with_headers(header,data,save_file);


    case 4 %Anxiety Word
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Word_After_Reward/' naming.subj_folder num2str(subj) '/LSSout_MNI_ANX.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/anx_after_reward/sub-' num2str(subj) '_REWARD_ANX.csv'];
        [patterns]=create_ROI_patterns(paths,mask_name,flags);

        %% Create Within Values
        disp('Creating Within Patterns For Anxiety Word Trials...')
        A=tril(corr(patterns),-1);%Within ROI, correlate all the exmemplars
        temp_r=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
        Within_R_Values=temp_r;

        %% Write Outputs
        header={'Within_Anxiety_Word_Trials'};
        data=Within_R_Values;
        save_file=[paths.save subj '_Within_Anxiety_Word_' mask_name '.txt'];
        save_data_with_headers(header,data,save_file);

    case 5 %Neutral Non-Reward
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Non_Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_non-reward_motion/sub-' num2str(subj) '.csv'];
        [patterns]=create_ROI_patterns(paths,mask_name,flags);

        %% Create Within Values
        disp('Creating Within Patterns For Neutral Non-Reward Trials...')
        A=tril(corr(patterns),-1);%Within ROI, correlate all the exmemplars
        temp_r=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
        Within_R_Values=temp_r;

        %% Write Outputs
        header={'Within_Neutral_Non_Reward_Trials'};
        data=Within_R_Values;
        save_file=[paths.save subj '_' mask_name '_Neutral_NoReward.txt'];
        save_data_with_headers(header,data,save_file);

    case 6 %Anxious Non-Reward
        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Non_Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_ANX.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/anxiety_non-reward_motion/sub-' num2str(subj) '.csv'];
        [patterns]=create_ROI_patterns(paths,mask_name,flags);

        %% Create Within Values
        disp('Creating Within Patterns For Anxious Non-Reward Trials...')
        A=tril(corr(patterns),-1);%Within ROI, correlate all the exmemplars
        temp_r=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
        Within_R_Values=temp_r;

        %% Write Outputs
        header={'Within_Anxious_Non_Reward_Trials'};
        data=Within_R_Values;
        save_file=[paths.save subj '_' mask_name '_Anxious_NoReward.txt'];
        save_data_with_headers(header,data,save_file);

    case 7 %Neutral & Anxiety Non-Reward

        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Non_Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_ANX.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/anxiety_non-reward_motion/sub-' num2str(subj) '.csv'];
        [anxiety_patterns]=create_ROI_patterns(paths,mask_name,flags);
        paths.data=[paths.top 'Data/Non_Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_non-reward_motion/sub-' num2str(subj) '.csv'];
        [neutral_patterns]=create_ROI_patterns(paths,mask_name,flags);

        %% Create Across Values
        disp('Creating Across Patterns...')

        anxiety_patterns=anxiety_patterns';
        neutral_patterns=neutral_patterns';
        A=size(anxiety_patterns,1);
        B=size(neutral_patterns,1);
        C=zeros(1,A*B);
        for a=1:A
            for b=1:B
                C(1,(a*b))=corr(anxiety_patterns(a,:)',neutral_patterns(b,:)');
            end
            progress(a,A,20);
        end
        temp_r=mean(atanh(C(find(C))));%Grab the non-zero values and do an r-to-z transform.
        Across_R_Values=temp_r;
        clear temp_r

        %% Write Outputs
        header={'Across_Non_Reward_Anxiety_and_Neutral_Trials'};
        data=Across_R_Values;
        save_file=[paths.save subj '_' mask_name '_Anxiety_vs_Neutral_NoReward.txt'];
        save_data_with_headers(header,data,save_file);

    case 8 %Neutral Reward & Neutral Non-Reward

        %% Get Patterns of Interest
        paths.data=[paths.top 'Data/Non_Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_non-reward_motion/sub-' num2str(subj) '.csv'];
        [non_reward_neutral_patterns]=create_ROI_patterns(paths,mask_name,flags);
        paths.data=[paths.top 'Data/Reward_Only_Trials/' naming.subj_folder num2str(subj) '/LSSout_MNI_NEU.nii'];
        paths.active_trials=[paths.reference '/Active_Trials/neutral_reward_motion/sub-' num2str(subj) '.csv'];
        [reward_neutral_patterns]=create_ROI_patterns(paths,mask_name,flags);

        %% Create Across Values
        disp('Creating Across Patterns...')

        non_reward_neutral_patterns=non_reward_neutral_patterns';
        reward_neutral_patterns=reward_neutral_patterns';
        A=size(non_reward_neutral_patterns,1);
        B=size(reward_neutral_patterns,1);
        C=zeros(1,A*B);
        for a=1:A
            for b=1:B
                C(1,(a*b))=corr(non_reward_neutral_patterns(a,:)',reward_neutral_patterns(b,:)');
            end
            progress(a,A,20);
        end
        temp_r=mean(atanh(C(find(C))));%Grab the non-zero values and do an r-to-z transform.
        Across_R_Values=temp_r;
        clear temp_r

        %% Write Outputs
        header={'Across_Non_Reward_Neutral_and_Reward_Neutral_Trials'};
        data=Across_R_Values;
        save_file=[paths.save subj '_' mask_name '_NonReward_Neutral_vs_Reward_Neutral.txt'];
        save_data_with_headers(header,data,save_file);

end

end