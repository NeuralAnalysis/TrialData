%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getMoveOnsetAndPeak(trial_data,min_ds)
%
%   This will find a time bin representing movement onset and peak speed
% Currently for kinematics but could be easily adapted to force, etc.
%
%   Constrains the looking to a window specified by start_idx and end_idx
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : parameter struct
%     .start_idx    : (string) field name of idx to use for start of window
%     .end_idx      : (string) field name to use for end of window
%     .which_method : (string) how to compute
%                           'peak' : uses acceleration and peak speed
%                           'thresh' : uses a basic velocity threshold
%                       Note: defaults to peak. In thresh, will not return
%                       a peak speed as a field.
%     .min_ds       : minimum diff(speed) to find movement onset
%     .s_thresh     : % speed threshold in cm/s (secondary method if first fails)
%     .peak_idx_offset  :  indices after start_idx to find max speed. 
%     .start_idx_offset :  indices after start_idx to find movement onset
%     .which_field     : which field to find movement onset from 
%     .field_idx       : idx of the above field to find movement onset from
%     .onset_name      : field to add for onset (default: 'movement_on')
%     .peak_name       : field to add for peak (default: 'peak_speed')
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
start_idx        =  'idx_go_cue';
end_idx          =  'idx_trial_end';
which_method     =  'peak';
min_ds           =  1.9;
s_thresh         =  10;
peak_idx_offset  = 0;
start_idx_offset = 0;
which_field      = 'speed';
field_idx        = 1;
onset_name       =  'movement_on';
peak_name        =  'peak_speed';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
verbose = false;
absolute_acc_thresh = [];
peak_divisor = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1, assignParams(who,params); end % overwrite defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end


for trial = 1:length(trial_data)
    % use velocity to find bin corresponding to movement onset, movement offset, and peak speed
    s = trial_data(trial).(which_field)(:,field_idx);
    
    % find the time bins where the monkey may be moving
    move_inds = false(size(s));
    move_inds(trial_data(trial).(start_idx)+start_idx_offset:trial_data(trial).(end_idx)) = true;
    
    [on_idx,peak_idx] = deal(NaN);
    if strcmpi(which_method,'peak')
        ds = [0; diff(s)];
        dds = [0; diff(ds)];
        peaks = [dds(1:end-1)>0 & dds(2:end)<0; 0];
        mvt_peak = find(peaks & (1:length(peaks))' > trial_data(trial).(start_idx)+peak_idx_offset & ds > min_ds & move_inds, 1, 'first');
        if ~isempty(mvt_peak)
            if isempty(absolute_acc_thresh)
                thresh = ds(mvt_peak)/peak_divisor; % Threshold is max of acceleration peak divided by divisor
            else
                thresh = absolute_acc_thresh;
            end
            on_idx = find(ds<thresh & (1:length(ds))'<mvt_peak & move_inds,1,'last');
            
            % check to make sure the numbers make sense
            if on_idx < trial_data(trial).(start_idx)+start_idx_offset
                % something is fishy. Fall back on threshold method
                warning('Something went wrong on trial %d...falling back to threshold method',trial)
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
            warning('Could not identify movement onset on trial %d',trial);
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
