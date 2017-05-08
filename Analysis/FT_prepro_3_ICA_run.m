%% Fieldtrip Preprocessing pipeline for RRD_EEG_2
% _____STEP 3.1: run ICA_____
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
sublist = [1:13 15:21];
dir_raw = 'Results/Raw/';
dir_ft = 'FTv3/';
path_in = [dir_ft 'prepro_2_1Hz/'];
file_in = 'OP1Hz_s';
condlist = {'REG10', 'REG10dev','RAND10', 'RAND10dev', 'dev_REG','dev_RAND'};   % condition labels
neeg = 128; nblocks = 6;
scount = 0;
% make some directories for saving intermediate data
mkdir([dir_ft 'ICA']);

%% Loop for automating ICA - returns structure called 'comp' with same attributes as 'data'
for s = sublist
    scount = scount + 1;
    load([path_in file_in num2str(s) '.mat']); % variable is Fieldtrp structure called 'data'

    %% ICA
    cfg        = [];
    nchanns = length(data.label);
    cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
%     cfg.runica.pca = nchanns - 1; % because rakn of matrix is reduced when referencing to average
    comp = ft_componentanalysis(cfg, data);
    writeFile  = [dir_ft 'ICA/ica' file_in num2str(s) '.mat'];
    save(writeFile, 'comp');
end
