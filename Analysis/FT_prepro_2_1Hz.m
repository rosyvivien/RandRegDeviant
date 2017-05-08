%% Fieldtrip Preprocessing pipeline for RRD_EEG_1
% _____STEP 2: Outlier Removal (1Hz ICA data)_____
% *** MUST be done after FT_prepro_2_sus & FT_prepro_2_dev because uses the
% this uses the outlier channels and trials defined in those steps
%
% _____PIPELINE_____
% FT_Prepro_1 -- read in, define trials, filters, downsample, append
% FT_Prepro_2 -- this file
% FT_Prepro_3 -- ICA (optional)
% FT_Prepro_4 -- Interpolate bad/missing channels, re-reference, average,
% LP filter at 30Hz
%
% _____STEPS_____

% - read in stored goodTrials & goodChanns and remove the rest (this
% includes removal of manual bad channels)
%____________ Rosy Southwell 2017-04____________


%% setup
close all; clearvars
sublist = [1];
dir_raw = 'Results/Raw/';
dir_ft = 'FTv3/';
path_in = [dir_ft 'prepro_1_1Hz/'];
file_in = 'P1Hz_s'; % stem of filename to read in (i.e. from last step of prepro 1
condlist = {'REG10', 'REG10dev','RAND10', 'RAND10dev', 'dev_REG','dev_RAND'};   % condition labels
neeg = 128;
nblocks = 6;
load badChannsManual.mat;
scount = 0;
% make some directories for saving intermediate data
mkdir([dir_ft 'prepro_2_1Hz']);
mkdir([dir_ft 'goodChanns']);

for s = sublist
    scount = scount + 1;
    load([path_in file_in num2str(s) '.mat']); % variable is Fieldtrp structure called 'data'
    
    load([dir_ft 'goodTrials_sus/sub' num2str(s)]);
    load([dir_ft 'goodChanns_sus/sub' num2str(s)]);
    goodChanns_sus = goodChanns(:);
    load([dir_ft 'goodChanns_dev/sub' num2str(s)]);
    goodChanns_dev = goodChanns(:);
    [tf,ix] = ismember(goodChanns_dev,goodChanns_sus);
    goodChanns = goodChanns_dev(tf);
    save([dir_ft 'goodChanns/sub' num2str(s)],'goodChanns');
    
    cfg = [];
    cfg.trials = goodTrials;
    cfg.channel = goodChanns;
    data = ft_selectdata(cfg,data);
    
    writeFile  = [dir_ft 'prepro_2_1Hz/' 'O' file_in num2str(s) '.mat'];
    save(writeFile, 'data');   
end