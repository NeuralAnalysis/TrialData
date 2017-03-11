%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function metric = getLearningMetrics(trial_data,params)
%
%   Computes and plot behavioral adaptation metrics.
%
% INPUTS:
%   trial_data: the struct
%   params: struct with the following options
%     .which_metric : which metric
%                       'angle': angular takeoff error
%                       'corr' : speed profile correlation coefficients
%                       'time' : time to target
%     .use_bl_ref   : (bool) whether to use diff from BL for angle
%     .time_window  : {'idx_start',bins after; 'idx_end',bins after} currently for angle
%     .corr_samples : how many datapoints to interpolate trajectory onto for corr
%     .vel_or_pos   : 'vel' or 'pos', for velocity or position (angle only)
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metric = getLearningMetrics(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
which_metric = 'angle';
time_window = {'idx_movement_on',0; 'idx_peak_speed',0};
use_bl_ref = true;
corr_samples = 1000;
vel_or_pos = 'vel';
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

utheta = unique([trial_data.target_direction]);
metric = zeros(length(trial_data),1);

switch lower(which_metric)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'angle'
        % get baseline error to each target
        bl_metric = zeros(length(utheta),1);
        if use_bl_ref
            for iDir = 1:length(utheta)
                bl_idx = getTDidx(trial_data,'epoch','bl','target_direction',utheta(iDir));
                if isempty(bl_idx)
                    error('No BL ref found for all targets');
                end
                
                temp_err = zeros(length(bl_idx),1);
                for iTrial = 1:length(bl_idx)
                    t1 = trial_data(bl_idx(iTrial)).(time_window{1,1})+time_window{1,2};
                    t2 = trial_data(bl_idx(iTrial)).(time_window{2,1})+time_window{2,2};
                    temp = trial_data(bl_idx(iTrial)).(vel_or_pos);
                    temp_err(iTrial) = angleDiff(minusPi2Pi(trial_data(bl_idx(iTrial)).target_direction), ...
                        atan2(temp(t2,2) - temp(t1,2), ...
                        temp(t2,1) - temp(t1,1)), ...
                        true,true);
                end
                bl_metric(iDir) = mean(temp_err);
            end, clear temp;
        end
        
        % get velocity at time of peak speed
        for iTrial = 1:length(trial_data)
            t1 = trial_data(iTrial).(time_window{1,1})+time_window{1,2};
            t2 = trial_data(iTrial).(time_window{2,1})+time_window{2,2};
            
            temp = trial_data(iTrial).(vel_or_pos);
            temp_err = angleDiff(minusPi2Pi(trial_data(iTrial).target_direction), ...
                atan2(temp(t2,2) - temp(t1,2), ...
                temp(t2,1) - temp(t1,1)), ...
                true,true);
            
            iDir = utheta==trial_data(iTrial).target_direction;
            metric(iTrial) = angleDiff(bl_metric(iDir),temp_err,true,true);
        end, clear temp;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'corr'
        % get baseline trace to each target
        bl_metric = zeros(length(utheta),corr_samples,2);
        for iDir = 1:length(utheta)
            bl_idx = getTDidx(trial_data,'epoch','bl','target_direction',utheta(iDir));
            bl_temp = zeros(length(bl_idx),2,corr_samples);
            if isempty(bl_idx)
                error('Corr needs a BL reference for each target');
            end
            
            for iTrial = 1:length(bl_idx)
                t1 = trial_data(bl_idx(iTrial)).(time_window{1,1})+time_window{1,2};
                t2 = trial_data(bl_idx(iTrial)).(time_window{2,1})+time_window{2,2};
                
                idx = t1:t2;
                temp = trial_data(bl_idx(iTrial)).vel;
                
                bl_temp(iTrial,1,:) = interp1(1:length(idx),temp(idx,1),linspace(1,length(idx),corr_samples));
                bl_temp(iTrial,2,:) = interp1(1:length(idx),temp(idx,2),linspace(1,length(idx),corr_samples));
            end
            bl_metric(iDir,:,:) = squeeze(mean(bl_temp,1))';
        end, clear temp bl_temp;
        
        for iTrial = 1:length(trial_data)
            t1 = trial_data(iTrial).(time_window{1,1})+time_window{1,2};
            t2 = trial_data(iTrial).(time_window{2,1})+time_window{2,2};
            
            idx = t1:t2;
            iDir = utheta==trial_data(iTrial).target_direction;
            
            temp = trial_data(iTrial).vel;
            
            temp = [interp1(1:length(idx),temp(idx,1),linspace(1,length(idx),corr_samples))', ...
                interp1(1:length(idx),temp(idx,2),linspace(1,length(idx),corr_samples))'];
            metric(iTrial) = corr2(squeeze(bl_metric(iDir,:,:)),temp)^2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'time'
        error('Time to target not implemented.');
        
    otherwise
        error('metric not recognized.');
end


