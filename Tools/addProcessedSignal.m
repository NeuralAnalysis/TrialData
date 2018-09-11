%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = addProcessedSignal(trial_data,params)
%
%   Adds spherical coordinates for opensim hand position to trial_data.
%   New signals are called:
%   'sphere_hand_pos'
%   'sphere_hand_vel'
%   'sphere_hand_acc'
%   Order of coordinates in new signals is: azimuth, elevation, r
%
% INPUTS:
%   trial_data : the struct
%   params : parameters struct
%       in_signals : input signals, specified in the format specified by
%           the output of check_signals
%       out_signals_name : what to name the output signals
%       processor : function handle to process in_signals to new signals
%
% OUTPUTS:
%   trial_data : the struct with .sphere_hand_* fields (reordered logically)
%
% Written by Raeed Chowdhury. Updated Jan 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = addProcessedSignal(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(isstruct(trial_data), 'First input must be trial_data struct!')
assert(isstruct(params), 'Second argument must be a paramters struct!')

in_signals = '';
out_signals_name = '';
processor = [];
assignParams(who,params)

% check params
in_signals = check_signals(trial_data,in_signals);
assert(ischar(out_signals_name) && ~isempty(out_signals_name),'Input a valid out_signals_name')
assert(isa(processor,'function_handle'),'processor has to be a function handle')

for trial = 1:length(trial_data)
    % loop over postfixes
    input_data = get_vars(trial_data(trial),in_signals);
    % process
    output_data = processor(input_data);

    % check output
    assert(size(input_data,1) == size(output_data,1),'Output data is different length from input, possible problem with processor?')

    % assign to new field
    trial_data(trial).(out_signals_name) = output_data;
end

% restore logical order
trial_data = reorderTDfields(trial_data);

