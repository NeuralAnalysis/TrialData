%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = addSphereHand2TD(trial_data)
%
%   Adds new polar-derived coordinates for body point position to trial_data.
%   New signals are called (for example)
%   'sph_hand_pos'
%   'sph_hand_vel'
%   'sph_hand_acc'
%   Order of coordinates in new signals for spherical is: azimuth, elevation, r
%   Order of coordinates in new signals for cylindrical is: theta, r, z
%
% INPUTS:
%   trial_data : the struct
%   params: param struct
%       method : either 'opensim' or 'markers'
%       point : either 'hand' or 'elbow' currently
%       coord : either 'sph' or 'cyl'
%
% OUTPUTS:
%   trial_data : the struct with .sphere_hand_* fields (reordered logically)
%
% Written by Raeed Chowdhury. Updated Jan 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = addCoordPoint2TD(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(isstruct(trial_data), 'First input must be trial_data struct!')
method = 'markers';
point = 'hand';
coord = 'sph';
if nargin>1
    assignParams(who,params);
end

% check inputs
assert(strcmpi(method,'markers') || strcmpi(method,'opensim'),'Method should be either ''markers'' or ''opensim''')
assert(strcmpi(point,'hand') || strcmpi(point,'elbow'),'Point should be either ''hand'' or ''elbow''')
assert(strcmpi(coord,'sph') || strcmpi(coord,'cyl'),'coord should be either ''sph'' or ''cyl''')

coord_postfix = {'pos','vel','acc'};
if strcmpi(method,'opensim')
    cart_prefix = {'X','Y','Z'};
    cart_postfix = strcat('_', point, {'Pos','Vel','Acc'});

    point_idx = cell(1,3);
    for postfix_ctr = 1:length(cart_postfix)
        [point_exists,point_idx{postfix_ctr}] = ismember(strcat(cart_prefix,cart_postfix(postfix_ctr)),trial_data(1).opensim_names);
        assert(all(point_exists),'Specified point %s does not exist in trial_data',point)
    end

    % duplicate signal name to match markers method organization
    cart_names = repmat({'opensim'},length(coord_postfix),1);

elseif strcmpi(method,'markers')
    % first get signal indices for point
    if strcmpi(point,'hand')
        markername = 'Marker_1';
    elseif strcmpi(point,'elbow')
        markername = 'Pronation_Pt1';
    else
        error('Specified point %s does not exist in trial_data',point)
    end
    [point_exists,point_idx] = ismember(strcat(markername,'_',{'x','y','z'}),trial_data(1).marker_names);
    assert(all(point_exists),'Specified point %s does not exist in trial_data',point)

    % duplicate point_idx to match opensim method organization
    point_idx = repmat({point_idx},length(coord_postfix),1);

    % fill out vel and acc fields if they don't exist yet
    if ~isfield(trial_data,'marker_vel')
        trial_data = smoothSignals(trial_data,struct('signals','markers'));
        trial_data = getDifferential(trial_data,struct('signals','markers','alias','marker_vel'));
    end
    if ~isfield(trial_data,'marker_acc')
        trial_data = smoothSignals(trial_data,struct('signals','marker_vel'));
        trial_data = getDifferential(trial_data,struct('signals','marker_vel','alias','marker_acc'));
    end

    cart_names = {'markers','marker_vel','marker_acc'};
else
    error('Unrecognized method to calculate spherical coordinates: %s',method)
end

% loop over postfixes to finally calculate
for postfix_ctr = 1:length(coord_postfix)
    % convert to new coordinate system
    if strcmpi(coord,'sph')
        trial_data = addProcessedSignal(trial_data,struct('in_signals',{{cart_names{postfix_ctr},point_idx{postfix_ctr}}},...
            'out_signals_name',['sph_hand_' coord_postfix{postfix_ctr}],...
            'processor',@cart2sph_wrap));
    elseif strcmpi(coord,'cyl')
        trial_data = addProcessedSignal(trial_data,struct('in_signals',{{cart_names{postfix_ctr},point_idx{postfix_ctr}}},...
            'out_signals_name',['cyl_hand_' coord_postfix{postfix_ctr}],...
            'processor',@cart2pol_wrap));
    else
        error('Unrecognized coordinate system')
    end
end

% restore logical order
trial_data = reorderTDfields(trial_data);

end

function sph = cart2sph_wrap(cart)
% helper function to collapse inputs and outputs
    [az, el, r] = cart2sph(cart(:,1),cart(:,2),cart(:,3));
    sph = [az el r];
end

function cyl = cart2pol_wrap(cart)
% helper function to collapse inputs and outputs
    [th,r,z] = cart2pol(cart(:,1),cart(:,2),cart(:,3));
    cyl = [th r z];
end
