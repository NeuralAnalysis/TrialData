%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function rw_td = getRWMovements(trial_data,params)
%
%   Separates random walk trials such that each movement in the trial is a
% unique trial_data entry. Can take some time before/after each go cue, but
% note that this will result in overlapping time points across trials. Will
% also compute a new movement onset and peak speed for each movement, if
% desired.
%
% Note: rewrites trial_id such that it is now [RW trial, movement number]
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .go_cue_name       : (string) name of go cue idx
%                               (default: 'idx_go_cue')
%       .extra_bins        : [TIME_BEFORE, TIME_AFTER] (in # bins)
%       .remove_incomplete : flag to remove trials where all targets were
%                               not acquired (e.g. nans in go cue)
%       .do_onset          : flag to compute movement onset (default true)
%
% OUTPUTS:
%   trial_data : struct separated by movements
%
% Written by Matt Perich. Updated March 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rw_td = getRWMovements(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS
go_cue_name        =  'idx_go_cue';
extra_bins         =  [0 0];
remove_incomplete  =  false;
do_onset           =  true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters to overwrite that aren't documented
start_name         =  'idx_trial_start';
end_name           =  'idx_trial_end';
onset_name         =  'idx_movement_on';
peak_name          =  'idx_peak_speed';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1, assignParams(who,params); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% take out trials where all targets weren't acquired
if remove_incomplete
    bad_trials = cellfun(@(x) ~any(isnan(x)),{trial_data.(go_cue_name)});
    trial_data(bad_trials) = [];
end

num_moves = sum(~isnan([trial_data.(go_cue_name)]));
rw_td = repmat(struct(),1,num_moves);

fn_time = getTDfields(trial_data,'time');
fn_meta = getTDfields(trial_data,'meta');

count = 0;
for trial = 1:length(trial_data)
    td = trial_data(trial);
    go_cues = td.(go_cue_name);
    
    for cue = 1:length(go_cues)
        if ~isnan(go_cues(cue))
            count = count + 1;
            % copy over the meta data
            for i = 1:length(fn_meta)
                rw_td(count).(fn_meta{i}) = td.(fn_meta{i});
            end
            % adjust trial_id
            rw_td(count).trial_id = [td.trial_id, cue];
            % target_center will be copied over in full, so index it
            rw_td(count).target_center = rw_td(count).target_center(cue,:);
            
            % copy over the time data
            
            idx_start = go_cues(cue)-extra_bins(1);
            if cue < length(go_cues)
                if isnan(go_cues(cue+1))
                    idx_end = td.(end_name)+extra_bins(2);
                else
                    idx_end = go_cues(cue+1)+extra_bins(2);
                end
            else
                idx_end = td.(end_name)+extra_bins(2);
            end
            
            % now add the new idx
            rw_td(count).(start_name) = extra_bins(1)+1;
            rw_td(count).(go_cue_name) = extra_bins(1)+1;
            rw_td(count).(end_name)   = idx_end - extra_bins(2) - go_cues(cue)+1;
            
            % check that the index won't crash
            if idx_end > size(td.(fn_time{1}),1)
                disp('Requested time extended beyond available trial data. Defaulting to last bin.');
                idx_end = size(td.(fn_time{1}),1);
            end
            % add time signals
            for i = 1:length(fn_time)
                temp = td.(fn_time{i});
                rw_td(count).(fn_time{i}) = temp(idx_start:idx_end,:);
            end
        end
    end
end

% get movement onset etc for each trial
if do_onset
    rw_td = getMoveOnsetAndPeak(rw_td,struct( ...
        'start_idx',start_name, ...
        'end_idx',end_name, ...
        'onset_name',onset_name, ...
        'peak_name',peak_name));
end

rw_td = reorderTDfields(rw_td);

