function [cfg] = ft_definetrial_chaitlab_RVS(cfg)
%Sort the trials in chronological order
% Definition of trials based on events at Ear Institute BIOSEMI system
% Trigger values are determined from the DURATION of the trigger pulses !
%
% FORMAT [trl, conditionlabels, S] = spm_eeg_definetrial(S)
% S                 - input structure (optional)
% (optional) fields of S:
%   S.D             - MEEG object or filename of M/EEG mat-file
%   S.timewin       - time window (in PST ms)
%   S.trialdef      - structure array for trial definition with fields (optional)
%       S.trialdef.conditionlabel - string label for the condition
%       S.trialdef.eventtype      - string (should be 'trigger_up')
%       S.trialdef.eventvalue     - trigger value for this condition
%       S.trialdef.trlshift       - shift the triggers by a fixed amount (sec)
%                                   (e.g. projector delay). One per
%                                   condition/trigegr type
%   S.filename      - path of BDF file to read
% OUTPUT:
%   trl             - Nx3 matrix [start end offset]
%   conditionlabels - Nx1 cell array of strings, label for each trial
%   S               - modified configuration structure (for history)
%__________________________________________________________________________
% Nicolas Barascud
% Adapted code from Vladimir Litvak, Robert Oostenveld
% Edited by Rosy Southwell - March 2017:
% - increased tolerance to +- 3
% - allow a different trialshift per condition
% - ignore trigger values / conditions which arent passed in in the trialdef, don't include
% these in cfg.event

%% Parameters
%--------------------------------------------------------------------------
% data = varargin{1};
% read the header, required to determine the stimulus channels and trial specification
hdr = ft_read_header(cfg.dataset);

begsample = 1;
endsample = hdr.nSamples*hdr.nTrials;
dataformat  = [];

% find the STATUS channel and read the values from it
schan = find(strcmpi(hdr.label,'STATUS'));
sdata = ft_read_data(cfg.dataset, 'header', hdr, 'dataformat', dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', schan);

pretrig  = cfg.trialdef.prestim;
posttrig = cfg.trialdef.poststim;

%% Read trigger channel (code adapted from ft_read_event.m)
%--------------------------------------------------------------------------

% trigger=sdata;

% convert to 32-bit integer representation and only preserve the lowest 24 bits
sdata = bitand(int32(sdata), 2^24-1);

byte1 = 2^8  - 1;
byte2 = 2^16 - 1 - byte1;
byte3 = 2^24 - 1 - byte1 - byte2;

% get the respective status and trigger bits
trigger = bitand(sdata, bitor(byte1, byte2)); % this is contained in the lower two bytes
trigger = bitand(sdata, byte1); % this is contained in the lower two bytes
epoch   = int8(bitget(sdata, 16+1));
cmrange = int8(bitget(sdata, 20+1));
battery = int8(bitget(sdata, 22+1));

% determine when the respective status bits go up or down
flank_trigger = diff([trigger]);
flank_epoch   = diff([0 epoch]);
flank_cmrange = diff([0 cmrange]);
flank_battery = diff([0 battery]);

%% Create event structure
%--------------------------------------------------------------------------

event       = [];
pad         = 0;
trigshift   = 0;
% convert the trigger into an event with a value at a specific sample
for i=find(flank_trigger>0)
    event(end+1).type   = 'trigger_up';        % distinguish between up and down flank
    event(end  ).sample = i + begsample;      % assign the sample at which the trigger has gone down
    event(end  ).value  = double(trigger(i+trigshift));      % assign the trigger value just _after_ going up
end

% Sort events in chronological order
[tmp, ind] = sort([event.sample]);
event = event(ind);

% Find distance between consecutive events
distance = find(flank_trigger<0) - find(flank_trigger>0);

% Remove events that are too close to each other
idx = find(distance< 4);
event(idx) = [];
distance(idx) = [];

% Find distance between consecutive events ? again
distance = double(distance);

figure(99); clf;
plot(trigger); hold on;
plot([event.sample], max(trigger),'o');
title('Trigger duration detection')
xlabel('Time (samples)')
ylabel('Trigger pulse amplitude')
ylim([double(min(trigger))-0 double(max(trigger))+0.1])
xlim([1 numel(trigger)])
legend({'Channel Amplitude' 'Detected Triggers'},'location','southoutside')

%% Safecheck. It's best to use MATLAB 2015a at this point.
if verLessThan('matlab','8.5.0')
    disp('Warning: No proper trigger safecheck was made.');
    disp('Matlab 2015a (or later) is required for optimal processing.');
    disp('Performing ''basic'' check instead (this is more likely to fail later...)');
    
    % Find unique trigger values in data (there should be as many as conditions)
    safecheck = unique(distance);
    
    if numel(safecheck) ~=numel(cfg.trialdef)
        disp([num2str(numel(safecheck)) ' unique trigger values were found in ' ]);
        disp(['trigger channel (' num2str(numel(cfg.trialdef)) ' required).']);
        disp('Attempting to continue anyway...');
    else
        disp(['A total ' num2str(numel(safecheck)) ' unique trigger values were found in ']);
        disp(['trigger channel (' num2str(numel(cfg.trialdef)) ' required). All OK !']);
    end
    
else % execute code for R2015a later
    % Find unique trigger values in data (there should be as many as conditions)
    tolerance = 2/max(distance); %% Tolerance decided
    safecheck = uniquetol(distance,tolerance);
    
    if numel(safecheck) ~=numel(cfg.trialdef.conditionlabel)
        disp([num2str(numel(safecheck)) ' unique trigger values were found in ']);
        disp(['trigger channel (' num2str(numel(cfg.trialdef.conditionlabel )) ' required). ']);
        %         return; %%Commented by Sijia: Retrive trials with only 2 out of 4
        %         triggers in EEG2 2016/02/09
    else
        disp(['A total ' num2str(numel(safecheck)) ' unique trigger values were found in']);
        disp(['trigger channel (' num2str(numel(cfg.trialdef.conditionlabel )) ' required). All OK !']);
    end
    
end

% Dirty fix for inaccurate trigger durations. Biosemi system codes event
% durations with ±1 time sample precision, so we manually correct these
% inaccuracies
tol = 3; % RVS changed from 2
for i=1:numel(cfg.trialdef.conditionlabel) % for all trial types
    % find trigger values ± 3 samples to be on the safe side
    idx = (distance >= (cfg.trialdef.eventvalue(i)-tol) ...
        & distance <= (cfg.trialdef.eventvalue(i)+tol)); %%%%% Tolerance
    distance(idx) = cfg.trialdef.eventvalue(i);
end

% Remove trigger_down events as we don't need them anymore
idx = strcmp({event.type},'trigger_down');
event(idx) = [];

for i=1:numel(event) % re-populate value with DURATION instead of AMPLITUDE of trigger event
    event(i).value = distance(i);
end
% remove any events not requested in trialdef
ix_rem = find(~ismember([event(:).value],cfg.trialdef.eventvalue));
if ~isempty(ix_rem)
event(ix_rem) = [];
distance(ix_rem) = [];
disp(['Removing ' num2str(length(ix_rem)) ' trials unrequested in trialdef; ' num2str(length(event)) ' remaining.']);
end
cfg.event = event;

%% Build trl matrix based on selected events
%--------------------------------------------------------------------------

for j=1:numel(cfg.trialdef.eventvalue)
    if ~isfield(cfg.trialdef,'trlshift')
        trlshift(j) = 0;
    elseif numel(cfg.trialdef.eventvalue) == numel(cfg.trialdef.trlshift)      
        trlshift(j) = round(cfg.trialdef.trlshift(j) * hdr.Fs); % assume passed as s
    elseif numel(cfg.trialdef.trlshift) == 1
        trlshift(j) = cfg.trialdef.trlshift * hdr.Fs;
    else
        error('cfg.trialdef.trlshift is wrong dimensions!')
    end
end
trl = [];
conditionlabels = {};
%% New
for i=1:numel(event)
    if ~strcmp(cfg.trialdef.eventtype,'trigger_up')
        disp('ERROR: S.trialdef.eventtype should be ''trigger_up''. Aborting!')
        return;
    end
    [icondition icondition]=find(cfg.trialdef.eventvalue==event(i).value);
    if isempty(icondition) % then we don't need this event as it is not in the list
    error('Wut. some unwanted trials are in the event structure...')    
    end
    trloff = round(pretrig*hdr.Fs); % assume passed as s
    trlbeg = event(i).sample - trloff; % trloff is prestim; positive number
    trldur = round((pretrig+posttrig)*hdr.Fs);% assume passed as s
    trlend = trlbeg + trldur;
    
    % Added by Rik in case wish to shift triggers (e.g, due to a delay
    % between trigger and visual/auditory stimulus reaching subject).
    trlbeg = trlbeg + trlshift(icondition);
    trlend = trlend + trlshift(icondition);
    
    % Add the beginsample, endsample and offset of this trial to the list
    trl = [trl; trlbeg trlend -trloff distance(i)]; % !! Fieldtrip difference prestim positive value means trial begins before trigger
    conditionlabels{end+1} = cfg.trialdef.conditionlabel(icondition);
end
cfg.trl = trl;
cfg.conditionlabels = conditionlabels;
