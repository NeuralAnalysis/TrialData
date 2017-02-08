function pds = getNeuronTuning(trial_data,which_method,params)
% will compute tuning curves using desired method
win = params.window; % {'idx_OF_START',BINS_AFTER; 'idx_OF_END', BINS_AFTER}
array = params.array;
dt = 0.01;

theta = [trial_data.target_direction]';

fr = zeros(length(trial_data),size(trial_data(1).([array '_spikes']),2));
for trial = 1:length(trial_data)
    temp = trial_data(trial).([array '_spikes']);
    idx = trial_data(trial).(win{1,1})+win{1,2}:trial_data(trial).(win{2,1})+win{2,2};
    fr(trial,:) = sum(temp(idx,:),1)./(dt*sum(idx));
end

tcs = regressTuningCurves(fr,theta,{'none'},'doparallel',false);

pds = tcs(:,3);