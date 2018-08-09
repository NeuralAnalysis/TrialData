%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = smoothSignals(trial_data,params)
%
% This function will smooth any signal.
%
% HINT: instead of params struct, can just pass SIGNALS input if you are
% okay just using the default parameters.
%
% INPUTS:
%   trial_data : the struct
%   params     :
%     .signals        : field names to smooth (single string or cell array)
%                          Note: must be passed in
%     .kernel_SD      : kernel s.d. for smoothing (default: 0.05)
%     .calc_rate      : flag to calculate rate (divide by bin size)
%                           Mostly meant for spiking data (default: false)
%     .sqrt_transform : flag to square root transform first (default: false)
%                           Mostly meant for spiking data
%
% OUTPUTS:
%   trial_data : the struct with all (signals{}) fields smoothed
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = smoothSignals(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
signals         =  []; 
kernel_SD       =  0.05;
calc_rate       =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters that aren't documented in the header
do_smoothing    =  true; % will just return trial_data if this is false
if nargin > 1
    if ~isstruct(params) % must be a signals input
        signals = params;
    else % overwrite defaults
        assignParams(who,params);
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end
if isempty(signals), error('Must provide one or more signals to smooth.'); end
if ~iscell(signals), signals = {signals}; end
bin_size = trial_data(1).bin_size;
% make sure the signal input formatting is good
signals = check_signals(trial_data(1),signals);
signals = signals(:,1); % don't need the idx if they exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if do_smoothing
    for trial = 1:length(trial_data)
        for i = 1:length(signals)
            data = trial_data(trial).(signals{i});
            if calc_rate, data = data./bin_size; end
            if do_smoothing
                data = smooth_data(data,bin_size,kernel_SD);
            end
            trial_data(trial).(signals{i}) = data;
        end
    end
end
