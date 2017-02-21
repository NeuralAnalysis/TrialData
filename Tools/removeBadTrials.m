%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data, bad_trials] = removeBadTrials(trial_data, params);
%
%   This function will identify and remove bad trials. Not much supported
% at the moment, but functionality can be expanded.
%       1) Will remove any trials with any idx_ fields that are NaN
%       2) Will remove trials according to ranges input
%
%
% INPUTS:
%   trial_data : trial data struct
%   params     : (struct) has parameter values
%       .ranges : {'idx_START','idx_END',[MIN_#BINS,MAX_#BINS]; etc...}
%                 ex: {'idx_go_cue','idx_movement_on',[5 30]} to remove
%                     reaction times smaller than 5 and larger than 30 bins
%
% OUTPUTS:
%   trial_data : struct with bad trials removed
%   bad_trials : struct containing said bad trials
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,bad_trials] = removeBadTrials(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ranges = [];
if nargin > 1, assignParams(who,params); end % overwrite defaults

bad_idx = false(1,length(trial_data));
for iTrial = 1:length(trial_data)
    err = false;
    
    td = trial_data(iTrial);
    
    % loop along all indices and make sure they aren't NaN
    fn = getTDfields(td,'idx');
    if any(cell2mat(cellfun(@(x) isnan(td.(x)),fn,'uni',0)))
        err = true;
    end
    
    %%%% LOOK FOR TRIALS THAT ARE OUTSIDE THE ALLOWABLE LENGTH
    if ~isempty(ranges)
        if size(ranges,2) ~= 3, error('Ranges input not properly formatted.'); end
        for i = 1:size(ranges,1)
            % define index values so I can check to make sure it's okay
            
            % If your requested values don't exist...
            if isempty(td.(ranges{i,1})) || isempty(td.(ranges{i,2}))
                error('idx references are outside trial range.');
            end
                
            if td.(ranges{i,2}) - td.(ranges{i,1}) < ranges{i,3}(1) || ...
                    td.(ranges{i,2}) - td.(ranges{i,1}) > ranges{i,3}(2)
                err = true;
            end
        end
    end
    
    if err, bad_idx(iTrial) = true; end
end
disp(['Pruning ' num2str(sum(bad_idx)) ' trials.']);
bad_trials = trial_data(bad_idx);
trial_data = trial_data(~bad_idx);