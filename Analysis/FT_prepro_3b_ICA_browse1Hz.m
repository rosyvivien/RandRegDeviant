%% Fieldtrip Preprocessing pipeline for RRD_EEG_2
% _____STEP 3.2: inspect ICA components & remove_____
% Run after outlier removal
% ***do this on the 1Hz high-passed data***
%
% _____PIPELINE_____
% FT_Prepro_1 -- read in, define trials, filters, downsample, append
% FT_Prepro_2 -- outlier removal
% FT_Prepro_3_ICArun -- this file
% FT_Prepro_4 -- Interpolate bad/missing channels, re-reference, average
%
% ____Rosy Southwell 2017-04________________

clearvars; close all
%% setup
close all; clearvars;
sublist = [12:13 15:21];
dir_raw = 'Results/Raw/';
dir_ft = 'FTv3/';
path_in = [dir_ft 'prepro_2_1Hz/'];
file_in = 'OP1Hz_s';
condlist = {'REG10', 'REG10dev','RAND10', 'RAND10dev', 'dev_REG','dev_RAND'};   % condition labels
neeg = 128; nblocks = 6;
load neighboursBiosemi128.mat;
scount = 0;
% make some directories for saving intermediate data
export = 1; plotdir = [dir_ft 'ICA/topo/']; mkdir(plotdir);

%% inspect & remove components - MANUAL STEP
ICAfile_in = ['ica' file_in ];
ICApath_in = [dir_ft 'ICA/' ];
for s = sublist 
    load([ICApath_in ICAfile_in num2str(s) '.mat']); % variable is Fieldtrp structure called 'comp'
    load([path_in file_in num2str(s) '.mat']); % variable is Fieldtrp structure called 'data'
    
    %% prepare layout
    cfg = [];
    cfg.layout = 'Biosemi128.lay';
    layout = ft_prepare_layout(cfg,data);
    
    %% Look in topo domain
    cfg = [];
    cfg.component = 1:20;       % specify the component(s) that should be plotted
    cfg.layout    = layout; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    f = figure(44); clf
    ft_topoplotIC(cfg, comp)
    set(f,'Position',[10 10 1700 1200]);
    suptitle(['Subject ' num2str(s)])
    if export
        set(f,'color',[1 1 1]);
        export_fig(f,'-dpng',[plotdir 'topo_' file_in  num2str(s) '.png'])
    end
    
    %% Look in time domain
    cfg = [];
    cfg.layout = layout; % specify the layout file that should be used for plotting
    cfg.viewmode = 'component';
    ft_databrowser(cfg, comp)
    title(['Subject ' num2str(s)])
    disp(['Subject ' num2str(s) '............'])
    %% place a breakpoint below.....
end
