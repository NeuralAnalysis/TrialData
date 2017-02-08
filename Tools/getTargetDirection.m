function trial_data = getTargetDirection(trial_data,params)
% compute time-varying target direction (relative to hand position) for
% each trials
reach_distance = 8; % cm
hold_time = 0.5; % s
dt = 0.01; % s
if nargin == 2
    if isfield(params,'reach_distance'), reach_distance = params.reach_distance; end
    if isfield(params,'dt'), dt = params.dt; end
    if isfield(params,'hold_time'), hold_time = params.hold_time; end
end

for i = 1:length(trial_data)
    tgt_pos = repmat([reach_distance*cos(trial_data(i).target_direction), reach_distance*sin(trial_data(i).target_direction)],size(trial_data(i).pos,1),1);
    pos = trial_data(i).pos;
    targ = atan2(tgt_pos(:,2)-pos(:,2),tgt_pos(:,1)-pos(:,1));
    
    % it's undefined once he enters outer target
    targ(trial_data(i).idx_trial_end - ceil(hold_time/dt):end) = 0;
    
    trial_data(i).targ = targ;
end