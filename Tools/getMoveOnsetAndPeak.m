%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getMoveOnsetAndPeak(trial_data,min_ds)
%
%   This will find a time bin representing movement onset and peak speed
% Currently for kinematics but could be easily adapted to force, etc.
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : parameter struct
%     .min_ds     : minimum diff(speed) to find movement onset
%     .s_thresh   : % speed threshold in cm/s (secondary method if first fails)
%
% OUTPUTS:
%   trial_data : same struct, with fields for:
%       idx_peak_speed
%       idx_movement_onset
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = getMoveOnsetAndPeak(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
min_ds     =  0.3;
s_thresh   =  7;
% these parameters aren't documented because I expect them to not need to
% change but you can overwrite them if you need to.
start_idx  =  'idx_go_cue';
end_idx    =  'idx_trial_end';
if nargin > 1, assignParams(who,params); end % overwrite defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for trial = 1:length(trial_data)
    % use velocity to find bin corresponding to movement onset, movement offset, and peak speed
    s = sqrt(trial_data(trial).vel(:,1).^2 + trial_data(trial).vel(:,2).^2);
    
    % peak will only be between go cue and trial end
    s = s(trial_data(trial).(start_idx):trial_data(trial).(end_idx));
    ds = [0; diff(s)];
    dds = [0; diff(ds)];
    peaks = [dds(1:end-1)>0 & dds(2:end)<0; 0];
    mvt_peak = find(peaks & (1:length(peaks))' > trial_data(trial).idx_go_cue & ds > min_ds, 1, 'first');
    
    if ~isempty(mvt_peak)
        thresh = ds(mvt_peak)/2;                             % Threshold is half max of acceleration peak
        on_idx = find(ds<thresh & (1:length(ds))'<mvt_peak,1,'last');
        % find movement peak as maximum velocity
        s(1:on_idx) = 0;
        [~, peak_idx] = max(s);
        
        % check to make sure the numbers make sense
        if on_idx <= trial_data(trial).idx_go_cue
            % something is fishy. Fall back on threshold method
            on_idx = find(s > s_thresh,1,'first');
        end
    else
        [on_idx,peak_idx] = deal(NaN);
    end
    
    trial_data(trial).idx_movement_on = on_idx;
    trial_data(trial).idx_peak_speed = peak_idx;
end

% restore logical order
trial_data = reorderTDfields(trial_data);

end