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
%   trial_data : the struct with differential field (reordered logically)
%
% TO DO:
%   - support for multiple signals?
%
% Written by Matt Perich. Updated Feb 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = getDifferential(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_value = 0;  % value for t = 0
signal   = ''; % signal to process
alias    = ''; % what to call the differentiated field
if nargin > 1
    assignParams(who,params);
else
    error('No parameters provided. Need to specify signal, at least.');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(alias)
    alias = ['d' signal];
end
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
    trial_data(trial).(alias) = [initial_value*ones(1,size(data,2)); diff(data)];
end

% restore logical order
trial_data = reorderTDfields(trial_data);

