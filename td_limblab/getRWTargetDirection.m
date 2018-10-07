%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getRWTargetDirection(trial_data, params)
% 
%   Computes angles for RW reaches for use for compatibility with VisData.
%   Takes in the structure, and finds the reach angle from one target to
%   the next. Note that there must be only single targets (run
%   getRWMovements first).
%
% INPUTS:
%   trial_data : the struct
%   params     : struct of parameters
%     
%
% OUTPUTS:
%   trial_data : struct with 'ReachDir' field added
% 
% Written by Matt Perich. Updated Feb 2017.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ trial_data ] = getRWTargetDirection(trial_data, params  )
start_idx_name = 'cursor_start';
target_pos_name = 'target_center';
if nargin > 1, assignParams(who,params); end % overwrite defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end

for i = 1:length(trial_data)
    pos = trial_data(i).(start_idx_name);
    tgt_pos = trial_data(i).(target_pos_name);
    targ = atan2(tgt_pos(:,2)-pos(:,2),tgt_pos(:,1)-pos(:,1));
        
    trial_data(i).target_direction = targ;
end

% restore logical order
trial_data = reorderTDfields(trial_data);

end

