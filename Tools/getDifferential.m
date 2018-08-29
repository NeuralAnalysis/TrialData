%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getDifferential(trial_data,signal,alias)
%
%   Differentiates a signal (e.g. get velocity from position)
%
% INPUTS:
%   trial_data : the struct
%   params     : struct input with the following fields
%       .signal       : (string) which signal (field) to use
%       .alias        : (string, optional) what name to use in struct
%                           Default is to add a 'd' in front of signal name
%       .intial_value : value for signal at t = 0
%
% OUTPUTS:
%   trial_data : the struct with differential field
%
% Written by Matt Perich. Updated Feb 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = getDifferential(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals   = ''; % signal to process
alias    = ''; % what to call the differentiated field
if nargin > 1
    if isstruct(params)
        assignParams(who,params);
    else % it's just the signal name
        signals = params;
    end
else
    error('No parameters provided. Need to specify signal, at least.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end
signals = check_signals(trial_data,signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(alias)
    alias = cell(size(signals,1),1);
    for iSig = 1:size(signals,1)
        alias{iSig} = ['d' signals{iSig,1}];
    end
elseif ischar(alias)
    alias = {alias};
else
    error('What kind of alias did you put in?')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(alias) ~= size(signals,1)
    error('Wrong number of aliases!');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSig = 1:size(signals,1)
    signal = signals{iSig,1};
    
    %   Check to see if field exists
    if ~isfield(trial_data,signal)
        error('Signal input is not a field of the trial_data struct.');
    end
    %   Check to ensure that requested signal is a time-varying field
    fn = getTDfields(trial_data,'time');
    if ~ismember(signal,fn)
        error('Signal input is not a time-varying field');
    end
    
    % loop along trials and differentiate
    for trial = 1:length(trial_data)
        data = trial_data(trial).(signal);
        [~,trial_data(trial).(alias{iSig})] = gradient(data,trial_data(trial).bin_size);
    end
end

% restore logical order
trial_data = reorderTDfields(trial_data);

