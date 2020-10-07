%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = filterSignals(trial_data,params)
%
% This function will filter any signal (using filtfilt).
%
% INPUTS:
%   trial_data : the struct
%   params     :
%     .signals        : field names to smooth (single string or cell array)
%                          Note: must be passed in
%     .filt_a          : A from filtfilt (denominator coefs of filter)
%     .filt_b          : B from filtfilt (numerator coefs of filter)
%
% OUTPUTS:
%   trial_data : the struct with all (signals{}) fields filtered
%
% Written by Raeed Chowdhury. Updated October 2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = filterSignals(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
signals         =  []; 
filt_a           =  1;
filt_b       =  1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters that aren't documented in the header
field_extra     =  '';   % if empty, defaults to input field name(s)
verbose = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if ~isstruct(params) % must be a signals input
        signals = params;
    else % overwrite defaults
        assignParams(who,params);
    end
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
  
for trial = 1:length(trial_data)
    for iSig = 1:size(signals,1)
        data = getSig(trial_data(trial),signals(iSig,1));
        data = filtfilt(filt_b,filt_a,data);
        trial_data(trial).([signals{iSig,1} field_extra{iSig}]) = data;
    end
end
