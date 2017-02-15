%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = smoothSpikes(trial_data,params)
%
% This function will smooth and/or square root transform spikes
%
% INPUTS:
%   trial_data : the struct
%   params     :
%     .sqrt_transform : flag to square root transform spikes (default: true)
%     .do_smoothing   : flag to convolve spikes with gaussian (default: true)
%     .kernel_SD      : kernel s.d. for smoothing (default: 2*bin_size)
%     .calc_fr        : flag to calculate FR (divide by bin size)
%
% OUTPUTS:
%   trial_data : the struct with all _spikes fields smoothed
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = smoothSpikes(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
sqrt_transform  =  true;
do_smoothing    =  true;
kernel_SD       =  0.05;
calc_fr         =  false;
if nargin > 1, assignParams(who,params); end % overwrite defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_size = trial_data(1).bin_size;

if do_smoothing || sqrt_transform % if you don't want to do either just passes back trial_data
    % find the _spikes fields
    fn_spikes = getTDfields(trial_data,'spikes');
    
    for trial = 1:length(trial_data)
        for i = 1:length(fn_spikes)
            fr = trial_data(trial).(fn_spikes{i});
            if sqrt_transform, fr = sqrt(fr); end
            if calc_fr, fr = fr./bin_size; end
            if do_smoothing
                fr = smoothSpikesForPCA(fr,bin_size,kernel_SD);
            end
            trial_data(trial).(fn_spikes{i}) = fr;
        end
    end
end