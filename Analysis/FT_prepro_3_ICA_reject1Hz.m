%% Fieldtrip Preprocessing pipeline for RRD_EEG_2
% _____STEP 3.3:  remove_____
% Run after outlier removal
% ***do this on the 1Hz high-passed data***
%
% _____PIPELINE_____
% FT_Prepro_1 -- read in, define trials, filters, downsample, append
% FT_Prepro_2 -- outlier removal
% FT_Prepro_3_ICA_reject1Hz -- this file
% FT_Prepro_4 -- Interpolate bad/missing channels, re-reference, average
%
% ____Rosy Southwell 2017-04________________

clearvars; close all
%% setup
close all; clearvars;
sublist = [1:13 15:21];
dir_raw = 'Results/Raw/';
dir_ft = 'FTv3/';
path_in = [dir_ft 'prepro_2_1Hz/'];
file_in = 'OP1Hz_s';
ICAfile_in = ['ica' file_in ];
ICApath_in = [dir_ft 'ICA/' ];

load neighboursBiosemi128.mat;
scount = 0;
% make some directories for saving intermediate data
mkdir([dir_ft 'prepro_3_1Hz']);
export = 1; plotdir = [dir_ft 'ICA/topo/']; mkdir(plotdir);

%% !! Please enter bad components using setBadComponents.m !!!
readyToRemove = 1;
setBadComponents;
load badComponents;
if readyToRemove
    for s = sublist % do in separate loop because this step requires human judgement!!
        load([path_in file_in num2str(s) '.mat']); % variable is Fieldtrp structure called 'data'
        load([ICApath_in ICAfile_in num2str(s) '.mat']); % variable is Fieldtrp structure called 'comp'
        
        % remove the bad components and backproject the data
        cfg = [];
        cfg.component = badComponents{s}; % to be removed component(s)
        data = ft_rejectcomponent(cfg, comp, data);
        
        %% save as I ...
        writeFile  = [dir_ft 'prepro_3_1Hz/' 'I' file_in num2str(s) '.mat'];
        save(writeFile, 'data');
    end
end