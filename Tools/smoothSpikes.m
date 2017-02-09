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
%     .bin_size       : size of time bins in trial_data (required for smoothing)
%     .kernel_SD      : kernel s.d. for smoothing (default: 2*bin_size)
%
% OUTPUTS:
%   trial_data : the struct with all _spikes fields smoothed
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = smoothSpikes(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(params,'sqrt_transform'), sqrt_transform = params.sqrt_transform; else, sqrt_transform = true; end
if isfield(params,'do_smoothing'), do_smoothing = params.do_smoothing; else, do_smoothing = true; end
if isfield(params,'bin_size'), bin_size = params.bin_size; else, bin_size = NaN; end
if isfield(params,'kernel_SD'), kernel_SD = params.kernel_SD; else, kernel_SD = 2*bin_size; end

if do_smoothing || sqrt_transform % if false just passes back trial_data
    % find the _spikes fields
    fn = fieldnames(trial_data(1));
    fn_spikes = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'_spikes')));
    
    for trial = 1:length(trial_data)
        for i = 1:length(fn_spikes)
            fr = trial_data(trial).(fn_spikes{i});
            if sqrt_transform, fr = sqrt(fr); end
            if do_smoothing
                fr = smoothSpikesForPCA(fr,bin_size,kernel_SD);
            end
            trial_data(trial).(fn_spikes{i}) = fr;
        end
    end
end