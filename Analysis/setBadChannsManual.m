% manual bad channels for RRD_EEG_2 as flagged during the experiment
% one entry for each subject
badChannsManual = cell(21,6);
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



