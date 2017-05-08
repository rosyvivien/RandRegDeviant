%% Fieldtrip Preprocessing pipeline for RRD_EEG_2
% _____STEP 3.3: use unmixing matrix from one dataset to remove bad components from separate
% dataset
% Run after outlier removal
% ***do this on the 0.1Hz high-passed data***
% read in the unmixing matrix from the ICA on 1Hz data
% estimate components timeseries from this
% reject bad components identified already from the 1Hz data
%
% _____PIPELINE_____
% FT_Prepro_1 -- read in, define trials, filters, downsample, append
% FT_Prepro_2 -- outlier removal
% FT_Prepro_3_doICA -- this file
% FT_Prepro_4 -- Interpolate bad/missing channels, re-reference, average
%
% ____Rosy Southwell 2017-04________________

clearvars; close all
%% setup
close all; clearvars;
sublist = [10:13 15:21];
dir_raw = 'Results/Raw/';
dir_ft = 'FTv3/';
path_in = [dir_ft 'prepro_2_sus/'];
file_in = 'OOP_s';
comp_in = 'OP1Hz_s'; % where to get unmixing matrix
ICAfile_in = ['ica' comp_in ];
ICApath_in = [dir_ft 'ICA/' ];
setBadComponents;
load badComponents;
load neighboursBiosemi128.mat;
scount = 0;
% make some directories for saving intermediate data
mkdir([dir_ft 'ICA']);
mkdir([dir_ft 'prepro_3_sus']);

for s = sublist % do in separate loop because this step requires human judgement!!
    load([ICApath_in ICAfile_in num2str(s) '.mat']); % variable is Fieldtrp structure called 'comp'
    load([path_in file_in num2str(s) '.mat']); % variable is Fieldtrp structure called 'data'
    
    %% remove extra bad channels identified in dev pipeline & save in the prepro_2 folder
    load([dir_ft 'goodChanns/sub' num2str(s)]);
    cfg = [];
    cfg.channel = goodChanns;
    data = ft_selectdata(cfg,data);
    
    writeFile  = [dir_ft 'prepro_2_sus/' 'O' file_in num2str(s) '.mat'];
    save(writeFile, 'data');
    
    %% estimate components using the unmixing matrix from ICA_run
    cfg=[];
    cfg.unmixing     = comp.unmixing;% NxN unmixing matrix
    cfg.topolabel    = comp.topolabel;%Nx1 cell-array with the channel labels
    comp = ft_componentanalysis(cfg, data);
    
    %% prepare layout
    cfg = [];
    cfg.layout = 'Biosemi128.lay';
    layout = ft_prepare_layout(cfg,data);
    
    %% Look in topo domain
    cfg = [];
    cfg.component = 1:20;       % specify the component(s) that should be plotted
    cfg.layout    = layout; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    cfg.markersize = 4;
    f = figure(44); clf
    ft_topoplotIC(cfg, comp)
    set(f,'Position',[10 10 1700 1200]);
    suptitle(['Subject ' num2str(s)])
    
    %% Look in time domain
    cfg = [];
    cfg.layout = layout; % specify the layout file that should be used for plotting
    cfg.viewmode = 'component';
    ft_databrowser(cfg, comp)
    title(['Subject ' num2str(s)])
    
    %% place a breakpoint below.....
    %% !! Please enter bad components using setBadComponents.m !!!
    readyToRemove = 1;
    
    if readyToRemove
        % remove the bad components and backproject the data
        cfg = [];
        cfg.component = badComponents{s}; % to be removed component(s)
        data = ft_rejectcomponent(cfg, comp, data);
        
        % save as I ...
        writeFile  = [dir_ft 'prepro_3_sus/' 'I' file_in num2str(s) '.mat'];
        save(writeFile, 'data');  
        clear data
        clear comp
    end
end

