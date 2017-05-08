%% Fieldtrip Preprocessing pipeline for RRD_EEG_1
% _____STEP 2: Outlier Removal (sustained response)_____
%
% _____PIPELINE_____
% FT_Prepro_1 -- read in, define trials, filters, downsample, append
% FT_Prepro_2 -- this file
% FT_Prepro_3 -- ICA (optional)
% FT_Prepro_4 -- Interpolate bad/missing channels, re-reference, average,
% LP filter at 30Hz
%
% _____STEPS_____
% EITHER
% - read in & remove badChannsManual
% - visual outlier removal (trials and channels)
% OR
% - read in stored goodTrials & goodChanns and remove the rest (this
% includes removal of manual bad channels)
%____________ Rosy Southwell 2017-04____________


%% setup
close all; clearvars
sublist = [3:13 15:21];
dir_raw = 'Results/Raw/';
dir_ft = 'FTv3/';
path_in = [dir_ft 'prepro_1_sus/'];
file_in = 'P_s'; % stem of filename to read in (i.e. from last step of prepro 1
condlist = {'REG10', 'REG10dev','RAND10', 'RAND10dev', 'dev_REG','dev_RAND'};   % condition labels
neeg = 128;
nblocks = 6;
load badChannsManual.mat;
scount = 0;
% make some directories for saving intermediate data
mkdir([dir_ft 'prepro_2_sus']);
mkdir([dir_ft 'goodTrials_sus']);
mkdir([dir_ft 'goodChanns_sus']);
rejectStoredOutliers = 0;

for s = sublist
    scount = scount + 1;
    load([path_in file_in num2str(s) '.mat']); % variable is Fieldtrp structure called 'data'
    
    switch rejectStoredOutliers
        case 0
            cfg = [];
            cfg.channel = 'EEG';
            for b = 1:6 % find channels bad in any blocks and remove from the visual outlier procedure
                bcm = badChannsManual{s,b};
                bc_all = {};
                bccount = 0;
                if ~isempty(bcm)
                    for c = 1:length(bcm)
                        bccount = bccount+1;
                        bc_all{bccount} = ['-' bcm{c}];
                    end
                end
            end
            okChanns = ft_channelselection({'EEG' bc_all{:}}, data.label) ;
            
            cfg.channel = okChanns;
            cfg.keepchannel = 'no'; % TODO where is the info stored about which channels rejected
            cfg.keeptrial = 'no';
            cfg.alim = [];
            cfg.latency = [-1 4];
            
            trl_old = data.trialinfo;
            
            data = ft_rejectvisual(cfg, data);
            goodTrials = data.trialinfo(:,2);
            goodChanns = data.label;
            save([dir_ft 'goodTrials_sus/sub' num2str(s)],'goodTrials');
            save([dir_ft 'goodChanns_sus/sub' num2str(s)],'goodChanns');
            
        case 1 % read in stored outliers & reject these
            load([dir_ft 'goodTrials_sus/sub' num2str(s)]);
            load([dir_ft 'goodChanns_sus/sub' num2str(s)]);
            cfg = [];
            cfg.trials = goodTrials;
            cfg.channel = goodChanns;
            data = ft_selectdata(cfg,data);
    end
    
    writeFile  = [dir_ft 'prepro_2_sus/' 'O' file_in num2str(s) '.mat'];
    save(writeFile, 'data');
    
end