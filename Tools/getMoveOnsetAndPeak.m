%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getMoveOnsetAndPeak(trial_data,min_ds)
%
%   This will find a time bin representing movement onset and peak speed
% Currently for kinematics but could be easily adapted to force
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   min_ds     : (optional) minimum diff(speed) to find movement onset
%
% OUTPUTS:
%   trial_data : same struct, with fields for:
%       idx_peak_speed and
%       idx_movement_onset
% 
% Written by Matt Perich. Updated Feb 2017.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = getMoveOnsetAndPeak(trial_data,min_ds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up a default value
if nargin < 2, min_ds = 0.3; end

for i = 1:length(trial_data)
    % use velocity to find bin corresponding to movement onset, movement offset, and peak speed
    s = sqrt(trial_data(i).vel(:,1).^2 + trial_data(i).vel(:,2).^2);
    
    % peak will only be between go cue and trial end
    s = s(trial_data(i).idx_go_cue:trial_data(i).idx_trial_end);
    ds = [0; diff(s)];
    dds = [0; diff(ds)];
    peaks = [dds(1:end-1)>0 & dds(2:end)<0; 0];
    mvt_peak = find(peaks & (1:length(peaks))' > trial_data(i).idx_go_cue & ds > min_ds, 1, 'first');
    
    if ~isempty(mvt_peak)
        thresh = ds(mvt_peak)/2;                             % Threshold is half max of acceleration peak
        on_idx = find(ds<thresh & (1:length(ds))'<mvt_peak,1,'last');
        % find movement peak as maximum velocity
        s(1:on_idx) = 0;
        [~, peak_idx] = max(s);
    else
        [on_idx,peak_idx] = deal(NaN);
    end
    trial_data(i).idx_movement_on = on_idx;
    trial_data(i).idx_peak_speed = peak_idx;
end

end