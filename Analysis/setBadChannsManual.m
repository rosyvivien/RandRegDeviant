% FIELDTRIP pipeline for RRD_EEG EEG experiments
% manual bad channels as flagged during the experiment
% Structure with a field for each subject and block.
% Each entry contains a cell array of strings which are the bad channels.

% _____EXAMPLES_____
% .single bad channel D8 for subject 1, block 2
% -- badChannsManual(1,2) = {{'D8'}};
% .multiple bad channels for subject 1, block 2
% -- badChannsManual(1,2) = {{'D8','A32'}};
% .multiple bad channels for subject 20, blocks 1:6
% -- badChannsManual(20,1:6) = repmat({{'A26','A27','A28'}},1,6);
% _________________________________________________________________
% Written by Rosy Southwell 2017-04

%% INITIALISE
badChannsManual = cell(21,6); % cell(nSubjects,nBlocks)

%% MANUALLY ENTER BAD CHANNELS
badChannsManual(2,2) = {{'D8'}};
badChannsManual(5,1) = {{'B20'}};
badChannsManual(5,4) = {{'C24'}};
badChannsManual(7,2) = {{'B3'}};
badChannsManual(7,6) = {{'B10'}};
badChannsManual(13,1:6) = repmat({{'B24'}},1,6);
badChannsManual(14,1:6) = repmat({{'B2'}},1,6);
badChannsManual(16,1:6) = repmat({{'B2'}},1,6);
badChannsManual(17,1:6) = repmat({{'A26', 'A27', 'A28', 'A29', 'A30', 'A31' ,'A32' ,'B2'}},1,6);
badChannsManual(19,1:6) = repmat({{'A28','A29','A30','A31','A32'}},1,6);
badChannsManual(20,1:6) = repmat({{'A26','A27','A28','A29','A30','A31','A32'}},1,6);
save('badChannsManual.mat','badChannsManual');

%% SAVE
save('badChannsManual.mat','badChannsManual');% manual bad channels for RRD_EEG_1 as flagged during the experiment


