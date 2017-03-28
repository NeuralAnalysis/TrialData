%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getMoveOnsetAndPeak(trial_data,min_ds)
%
%   This will find a time bin representing movement onset and peak speed
% Currently for kinematics but could be easily adapted to force, etc.
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : parameter struct
%     .which_method : (string) how to compute
%                           'peak' : uses acceleration and peak speed
%                           'thresh' : uses a basic velocity threshold
%                       Note: defaults to peak. In thresh, will not return
%                       a peak speed as a field.
%     .min_ds       : minimum diff(speed) to find movement onset
%     .s_thresh     : % speed threshold in cm/s (secondary method if first fails)
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
which_method  =  'peak';
min_ds        =  0.3;
s_thresh      =  7;
% these parameters aren't documented because I expect them to not need to
% change but you can overwrite them if you need to.
start_idx     =  'idx_go_cue';
end_idx       =  'idx_trial_end';
onset_name    =  'movement_on';
peak_name     =  'peak_speed';
if nargin > 1, assignParams(who,params); end % overwrite defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some pre-processing
td = getSpeed(trial_data);

for trial = 1:length(trial_data)
    % use velocity to find bin corresponding to movement onset, movement offset, and peak speed
    s = td(trial).speed;
    
    % find the time bins where the monkey may be moving
    move_inds = false(size(s));
    move_inds(td(trial).(start_idx):td(trial).(end_idx)) = true;
    
    [on_idx,peak_idx] = deal(NaN);
    if strcmpi(which_method,'peak')
        ds = [0; diff(s)];
        dds = [0; diff(ds)];
        peaks = [dds(1:end-1)>0 & dds(2:end)<0; 0];
        mvt_peak = find(peaks & (1:length(peaks))' > td(trial).(start_idx) & ds > min_ds & move_inds, 1, 'first');
        if ~isempty(mvt_peak)
            thresh = ds(mvt_peak)/2; % Threshold is half max of acceleration peak
            on_idx = find(ds<thresh & (1:length(ds))'<mvt_peak & move_inds,1,'last');
            
            % check to make sure the numbers make sense
            if on_idx <= td(trial).(start_idx)
                % something is fishy. Fall back on threshold method
                on_idx = NaN;
            end
        end
        % peak is max velocity during movement
        temp = s; temp(~move_inds) = 0;
        [~,peak_idx] = max(temp);
    end
    
    if isempty(on_idx) || isnan(on_idx)
        on_idx = find(s > s_thresh & move_inds,1,'first');
        if isempty(on_idx) % usually means it never crosses threshold
            warning('Could not identify movement onset');
            on_idx = NaN;
        end
    end
    trial_data(trial).(['idx_' onset_name]) = on_idx;
    if strcmpi(which_method,'peak')
        trial_data(trial).(['idx_' peak_name]) = peak_idx;
    end
end

% restore logical order
trial_data = reorderTDfields(trial_data);

end