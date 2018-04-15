%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = addCylHand2TD(trial_data)
%
%   Adds cylindrical coordinates for opensim hand position to trial_data.
%   New signals are called:
%   'cyl_hand_pos'
%   'cyl_hand_vel'
%   'cyl_hand_acc'
%   Order of coordinates in new signals is: azimuth, elevation, r
%
% INPUTS:
%   trial_data : the struct
%
% OUTPUTS:
%   trial_data : the struct with .cyl_hand_* fields (reordered logically)
%
% Written by Raeed Chowdhury. Updated Jan 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = addCylHand2TD(trial_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cart_prefix = {'X','Y','Z'};
cart_postfix = {'_handPos','_handVel','_handAcc'};
cyl_postfix = {'pos','vel','acc'};
for trial = 1:length(trial_data)
    % check for opensim field
    if isfield(trial_data(trial),'opensim') && isfield(trial_data(trial),'opensim_names')
        % loop over postfixes
        for postfix_ctr = 1:length(cart_postfix)
            % check for handPos
            if sum(contains(trial_data(trial).opensim_names,cart_postfix{postfix_ctr})) == 3
                % convert to spherical
                coords = [];
                for prefix_ctr = 1:length(cart_prefix)
                    coords = [coords trial_data(trial).opensim(:,contains(trial_data(trial).opensim_names,...
                                                                        [cart_prefix{prefix_ctr} cart_postfix{postfix_ctr}]))];
                end
                [az,el,r] = cart2pol(coords(:,1),coords(:,2),coords(:,3));
                trial_data(trial).(['cyl_hand_' cyl_postfix{postfix_ctr}]) = [az el r];
            else
                warning(['no hand marker ' cyl_postfix{postfix_ctr} ' found for one or more trials'])
            end
        end
    else
        error('No opensim information in one or more trials')
    end
end

% restore logical order
trial_data = reorderTDfields(trial_data);

