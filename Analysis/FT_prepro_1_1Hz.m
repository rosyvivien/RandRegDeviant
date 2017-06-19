%% Fieldtrip Preprocessing pipeline for RRD_EEG_1
%% 1HZ for ICA
% Run before outliers, jumps, ICA, interp & rereference, average
%
% _____PIPELINE_____
% FT_Prepro_1 -- this file
% FT_Prepro_2 -- outliers & manual bad channel removal
% FT_Prepro_3 -- ICA (optional)
% FT_Prepro_4 -- Interpolate bad/missing channels, re-reference, average%
%
% _____STEPS_____
% HP 1 Hz
% define & apply trial
% LP 100 Hz
% downsample 200
% merge
% bc
%____________ Rosy Southwell 2017-03____________

%% setup
clearvars; close all
sublist = [1:13 15:21];
dir_raw = 'Results/Raw/';
dir_ft = 'FTv3/'; mkdir(dir_ft); %save intermediate files here
condlist = {'REG10', 'REG10dev','RAND10', 'RAND10dev', 'dev_REG','dev_RAND'};   % condition labels
triglist    = [50 60 70 80 100 120]; % list of triggers (in the same order as conditions)
nblocks = 6;
scount = 0;
% make some directories for saving intermediate data
mkdir([dir_ft 'prepro_1_1Hz/']);

for s = sublist
    scount = scount +1;
    disp(['Subject ' num2str(scount) ' of ' num2str(length(sublist))])
    for b = 1:nblocks
        file_raw = [dir_raw 'sub' num2str(s) '/sub' num2str(s) '_b' num2str(b) '.bdf'];
        
        %% Define 'dummy' trial which is all data: needed before filtering
        %         but we don't want to use actual epochs as our cutoff frequency is
        %         very low (0.1Hz) - thus very long data segments needed to avoid edge artefacts.
        cfg = [];
        cfg.dataset =  file_raw;
        cfg.trialdef.triallength = Inf;
        cfg.trialdef.ntrials = 1;
        cfg = ft_definetrial(cfg);
        
        %% High Pass 0.1Hz
        HPf = 0.1;
        cfg.demean = 'yes';
        cfg.hpfilter = 'yes';
        cfg.hpfreq = HPf;
        cfg.hpfiltord = 3;
        longdata = ft_preprocessing(cfg);
        
        %% Define trial (has to be done before ft_preprocessing
        cfg = [];
        cfg.dataset =  file_raw; % your filename with file extension;
        cfg.trialdef.eventtype  = 'trigger_up'; % Status notation maybe Biosemi's  tigger name
        cfg.trialdef.eventvalue = triglist; % your event value
        cfg.trialdef.conditionlabel = condlist;
        cfg.trialdef.prestim    = 1;  % before stimulation (sec), positive value means trial begins before trigger
        cfg.trialdef.poststim   = 4; % after stimulation (sec) , only use positive value
        [cfg] = ft_definetrial_chaitlab_RVS(cfg);
        
        if length(cfg.event)>150 %in case accidentally left recording when restarted block
            cfg.event = cfg.event(end-150+1:end);
            cfg.trl = cfg.trl(end-150+1:end,:);
            cfg.event = cfg.conditionlabels(end-150+1:end);
        end
        
        % apply new epochs
        data = ft_redefinetrial(cfg, longdata);
        
        %        Find possible NaNs in the data
        findNaNs
        cfg = [];
        cfg.trials = keepTrials;
        
        %% Low pass 100 % required before downsample & needs to be at most 0.5* downsamplefs because of aliasing
        LPf = 100;
        cfg.lpfilter = 'yes';
        cfg.lpfreq = LPf;
        cfg.lpfiltord = 5;
        data = ft_preprocessing(cfg,data);
        
        %% Downsample
        resamplefs = 200;
        cfg.resamplefs = resamplefs;
        cfg.detrend = 'no';
        data = ft_resampledata(cfg,data);
        blockData{b} = data;
    end
    
    %% Merge
    cfg = [];
    data = ft_appenddata(cfg,blockData{1}, blockData{2}, blockData{3}, blockData{4}, blockData{5}, blockData{6});
    data.trialinfo(:,2) = 1:length(data.trialinfo); % append a column of info for trial number
    blockData={};
    
    %% Baseline Correct
    cfg = [];
    cfg.demean = 'yes';
    cfg.baselinewindow = [-1 0];
    data = ft_preprocessing(cfg, data);
    writeFile  = [dir_ft 'prepro_1_1Hz/P1Hz_s' num2str(s) '.mat'];
    save( writeFile, 'data');
end