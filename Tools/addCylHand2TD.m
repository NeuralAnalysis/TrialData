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
function trial_data = addCylHand2TD(trial_data,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(isstruct(trial_data), 'First input must be trial_data struct!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
    method = 'marker';
end

if strcmpi(method,'opensim')
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
                    % convert to cylindrical
                    coords = [];
                    for prefix_ctr = 1:length(cart_prefix)
                        coords = [coords trial_data(trial).opensim(:,contains(trial_data(trial).opensim_names,...
                                                                            [cart_prefix{prefix_ctr} cart_postfix{postfix_ctr}]))];
                    end
                    [th,r,z] = cart2pol(coords(:,1),coords(:,2),coords(:,3));
                    trial_data(trial).(['cyl_hand_' cyl_postfix{postfix_ctr}]) = [th r z];
                else
                    warning(['no hand marker ' cyl_postfix{postfix_ctr} ' found for one or more trials'])
                end
            end
        else
            error('No opensim information in one or more trials')
        end
    end
elseif strcmpi(method,'markers')
    hand_idx = 1:3;
    cyl_postfix = {'pos','vel','acc'};
    if ~isfield(trial_data,'marker_vel')
        trial_data = smoothSignals(trial_data,struct('signals','markers'));
        trial_data = getDifferential(trial_data,struct('signals','markers','alias','marker_vel'));
    end
    if ~isfield(trial_data,'marker_acc')
        trial_data = smoothSignals(trial_data,struct('signals','marker_vel'));
        trial_data = getDifferential(trial_data,struct('signals','marker_vel','alias','marker_acc'));
    end
    cart_names = {'markers','marker_vel','marker_acc'};
    for trial = 1:length(trial_data)
        % loop over postfixes
        for postfix_ctr = 1:length(cyl_postfix)
            % check for opensim field
            if isfield(trial_data(trial),cart_names{postfix_ctr})
                % convert to cylindrical
                coords = get_vars(trial_data(trial),{cart_names{postfix_ctr},hand_idx});
                [th,r,z] = cart2pol(coords(:,1),coords(:,2),coords(:,3));
                trial_data(trial).(['cyl_hand_' cyl_postfix{postfix_ctr}]) = [th r z];
            else
                error('No opensim information in one or more trials')
            end
        end
    end
else
    error('Unrecognized method to calculate cylindrical coordinates: %s',method)
end

% restore logical order
trial_data = reorderTDfields(trial_data);
