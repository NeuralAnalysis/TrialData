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
%     .width          : kernel width (std. dev) for smoothing (default: 0.05)
%     .calc_rate      : flag to calculate rate (divide by bin size)
%                           Mostly meant for spiking data (default: false)
%
% OUTPUTS:
%   trial_data : the struct with all (signals{}) fields smoothed
%
% Written by Matt Perich. Updated March 2019.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = smoothSignals(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
signals         =  []; 
width           =  0.05;
calc_rate       =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters that aren't documented in the header
field_extra     =  '';   % if empty, defaults to input field name(s)
do_smoothing    =  true; % will just return trial_data if this is false
verbose = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if ~isstruct(params) % must be a signals input
        signals = params;
    else % overwrite defaults
        assignParams(who,params);
    end
end 

if isfield(params,'kernel_SD')
    warning('smoothSignals changed!!! the kernel_SD input is now called width. it looks like you passed in kernel_SD. Change it to width.');
    width = params.kernel_SD;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
bin_size = trial_data(1).bin_size;
% make sure the signal input formatting is good
signals = check_signals(trial_data(1),signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check output field addition
field_extra  = check_field_extra(field_extra,signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
if do_smoothing
    for trial = 1:length(trial_data)
        for iSig = 1:size(signals,1)
            data = getSig(trial_data(trial),signals(iSig,:));
            if calc_rate, data = data./bin_size; end
            if do_smoothing
                data = smooth_data(data,bin_size,width);
            end
            trial_data(trial).([signals{iSig,1} field_extra{iSig}]) = data;
        end
    end
end
