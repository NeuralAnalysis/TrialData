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
%     .target_dir_fieldname : name for target direction field (default: target_direction)
%
% Written by Matt Perich. Updated Jan 2020 by Raeed Chowdhury.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metric = getLearningMetrics(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
which_metric = 'angle';
time_window = {'idx_movement_on',0; 'idx_peak_speed',0};
use_bl_ref = true;
fit_bl_ref_curve = false;
corr_samples = 1000;
vel_or_pos = 'vel';
target_dir_fieldname = 'target_direction';
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end

utheta = unique([trial_data.(target_dir_fieldname)]);
metric = zeros(length(trial_data),1);

switch lower(which_metric)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'angle'
        % get baseline error to each target
        bl_metric = zeros(length(utheta),1);
        if use_bl_ref
            if fit_bl_ref_curve
                % try fitting sins and cosines for bl_ref
                num_freqs = 3;
                [bl_idx,bl_trials] = getTDidx(trial_data,'epoch','bl');
                theta = [bl_trials.(target_dir_fieldname)];
                sincos = [sin(theta'*(1:num_freqs)) cos(theta'*(1:num_freqs))];
                usincos = [sin(utheta'*(1:num_freqs)) cos(utheta'*(1:num_freqs))];
                if isempty(bl_idx)
                    error('No BL ref found');
                end
                
                temp_err = zeros(length(bl_idx),1);
                for iTrial = 1:length(bl_idx)
                    t1 = trial_data(bl_idx(iTrial)).(time_window{1,1})+time_window{1,2};
                    t2 = trial_data(bl_idx(iTrial)).(time_window{2,1})+time_window{2,2};
                    temp = trial_data(bl_idx(iTrial)).(vel_or_pos);
                    temp_err(iTrial) = angleDiff(minusPi2Pi(trial_data(bl_idx(iTrial)).(target_dir_fieldname)), ...
                        atan2(temp(t2,2) - temp(t1,2), temp(t2,1) - temp(t1,1)), ...
                        true,true);
                end

                % fit function
                lm_bl = fitlm(sincos,temp_err);
                bl_metric = lm_bl.predict(usincos);
            % end
            % if true
            else
                for iDir = 1:length(utheta)
                    bl_idx = getTDidx(trial_data,'epoch','bl',target_dir_fieldname,utheta(iDir));
                    if isempty(bl_idx)
                        error('No BL ref found for all targets');
                    end
                    
                    temp_err = zeros(length(bl_idx),1);
                    for iTrial = 1:length(bl_idx)
                        t1 = trial_data(bl_idx(iTrial)).(time_window{1,1})+time_window{1,2};
                        t2 = trial_data(bl_idx(iTrial)).(time_window{2,1})+time_window{2,2};
                        t2 = min(t2,trial_data(bl_idx(iTrial)).idx_endTime);
                        
                        temp = trial_data(bl_idx(iTrial)).(vel_or_pos);
                        temp_err(iTrial) = angleDiff(minusPi2Pi(trial_data(bl_idx(iTrial)).(target_dir_fieldname)), ...
                            atan2(temp(t2,2) - temp(t1,2), temp(t2,1) - temp(t1,1)), ...
                            true,true);
                    end
                    bl_metric(iDir) = mean(temp_err);
                end, clear temp;
                % figure; plot(utheta,bl_metric','o');
                % theta = linspace(0,2*pi,100);
                % sincos = [sin(theta'*(1:num_freqs)) cos(theta'*(1:num_freqs))];
                % blspline = lm_bl.predict(sincos);
                % hold on
                % plot(theta,blspline)
            end
        end
        
        % get velocity at time of peak speed
        for iTrial = 1:length(trial_data)
            t1 = trial_data(iTrial).(time_window{1,1})+time_window{1,2};
            t2 = trial_data(iTrial).(time_window{2,1})+time_window{2,2};
            t2 = min(t2,trial_data(iTrial).idx_endTime);
            
            temp = trial_data(iTrial).(vel_or_pos);
            temp_err = angleDiff(minusPi2Pi(trial_data(iTrial).(target_dir_fieldname)), ...
                atan2(temp(t2,2) - temp(t1,2), temp(t2,1) - temp(t1,1)), ...
                true,true);
            
            iDir = utheta==trial_data(iTrial).(target_dir_fieldname);
            metric(iTrial) = angleDiff(bl_metric(iDir),temp_err,true,true);
        end, clear temp;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'corr'
        % get baseline trace to each target
        bl_metric = zeros(length(utheta),corr_samples,2);
        for iDir = 1:length(utheta)
            bl_idx = getTDidx(trial_data,'epoch','bl',target_dir_fieldname,utheta(iDir));
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
            iDir = utheta==trial_data(iTrial).(target_dir_fieldname);
            
            temp = trial_data(iTrial).vel;
            
            temp = [interp1(1:length(idx),temp(idx,1),linspace(1,length(idx),corr_samples))', ...
                interp1(1:length(idx),temp(idx,2),linspace(1,length(idx),corr_samples))'];
            metric(iTrial) = corr2(squeeze(bl_metric(iDir,:,:)),temp)^2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'curvature'
        % get baseline curvature to each target
        bl_metric = zeros(length(utheta),1);
        if use_bl_ref
            for iDir = 1:length(utheta)
                bl_idx = getTDidx(trial_data,'epoch','bl',target_dir_fieldname,utheta(iDir));
                if isempty(bl_idx)
                    error('No BL ref found for all targets');
                end
                
                temp_err = zeros(length(bl_idx),1);
                for iTrial = 1:length(bl_idx)
                    t1 = trial_data(bl_idx(iTrial)).(time_window{1,1})+time_window{1,2};
                    t2 = trial_data(bl_idx(iTrial)).(time_window{2,1})+time_window{2,2};
                    temp = trial_data(bl_idx(iTrial)).(vel_or_pos);
                    temp_err(iTrial) = median(curvature(temp(t1:t2,:)));
                end
                bl_metric(iDir) = mean(temp_err);
            end, clear temp;
        end
        
        % 
        for iTrial = 1:length(trial_data)
            t1 = trial_data(iTrial).(time_window{1,1})+time_window{1,2};
            t2 = trial_data(iTrial).(time_window{2,1})+time_window{2,2};
            
            temp = trial_data(iTrial).(vel_or_pos);
            temp_err = median(curvature(temp(t1:t2,:)));
            
            iDir = utheta==trial_data(iTrial).(target_dir_fieldname);
            metric(iTrial) = temp_err - bl_metric(iDir);
        end, clear temp;
        
        metric = -metric;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'time'
        error('Time to target not implemented.');
        
    otherwise
        error('metric not recognized.');
end


