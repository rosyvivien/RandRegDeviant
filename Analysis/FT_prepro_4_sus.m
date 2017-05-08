%% Fieldtrip Preprocessing pipeline for RRD_EEG_1
%
% _____PIPELINE_____
% FT_Prepro_1 --
% FT_Prepro_2 -- outliers & jumps
% FT_Prepro_3 -- ICA (optional)
% FT_Prepro_4 -- Interpolate bad/missing channels, re-reference, baseline, average
%
%____________ Rosy Southwell 2017-03____________
clearvars; close all
%% setup
sublist = [1:13 15:21];
dir_raw = 'Results/Raw/';
dir_ft = 'FTv3/';
condlist = {'REG10', 'REG10dev','RAND10', 'RAND10dev'};   % condition labels
triglist    = [50 60 70 80]; % list of triggers (in the same order as conditions)
load badChannsManual.mat;
load neighboursBiosemi128.mat;
scount = 0;
% make some directories for saving intermediate data
mkdir([dir_ft 'prepro_4_sus']);
mkdir([dir_ft 'average_sus/']);
path_in = [dir_ft 'prepro_3_sus/'];
LP = 1; % apply low-pass (before averaging)??

for s = sublist
    file_in = 'IOP_s';
    load([path_in file_in num2str(s) '.mat']); % variable is Fieldtrp structure called 'data'
    
    %% Interpolate over bad channels
    cfg = [];
    cfg.layout = 'biosemi128.lay';
    cfg.neighbours = neighbours;
    % which channels are missing?
    allEEG = {neighbours.label};
    missing = setdiff(allEEG,data.label);
    cfg.missingchannel = missing;
    
    data = ft_channelrepair(cfg,data);
    if length(data.label)<128
        disp(['nChanns ' num2str(length(data.label))])
        pause
    end
    
    %% Re-reference
    cfg.reref = 'yes';
    cfg.refchannel =  'all';
    if isfield(cfg,'resamplefs')
        cfg = rmfield(cfg, 'resamplefs');
    end
    data = ft_preprocessing(cfg,data);
    %     writeFile  = [dir_ft 'Prepro4/C' file_in num2str(s) '.mat'];
    %     save(writeFile, 'data');
    
    %% BASELINE CORRECT
    cfg = [];
    cfg.demean = 'yes';
    cfg.baselinewindow = [-1 0];
    data = ft_preprocessing(cfg, data);
    writeFile  = [dir_ft 'prepro_4_sus/BC' file_in num2str(s) '.mat'];
    save( writeFile, 'data');
    
    %% average over trials
    data_all = data; % put aside separate variable for all conditions
    for c = 1:length(condlist)
        cfg = [];
        cfg.trials = find(data_all.trialinfo(:,1) == triglist(c));
        cfg.removemean='no';
        cfg.keeptrials = 'no';
        data = ft_timelockanalysis(cfg,data_all);
        %         writeFile  = [dir_ft 'AveLite/' condlist{c} '_ABC' file_in num2str(s) '.mat'];
        %         save(writeFile, 'data');
        if LP
            %% LP filter
            cfg = [];
            LPf = 30;
            cfg.lpfilter = 'yes';
            cfg.lpfreq = LPf;
            cfg.lpfiltord = 5;
            data = ft_preprocessing(cfg,data);
            
        end
        eval(['aveData.' condlist{c} ' = data']);
        
    end
    writeFile  = [dir_ft 'average_sus/' 'F30ABC' file_in num2str(s) '.mat'];
    save( writeFile, 'aveData');
end