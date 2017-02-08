%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot GPFA results
%   vis_data is focused on single trial plotting, including GPFA
%   trajectories, whereas this function will be intended for summarizing
%   across trials
%
% Can also add kinematic traces grouped across trials. Can't add any
%   spiking activity right now because grouping across trials doesn't make
%   sense unless we are smoothing it into a continuous firing rate
%
% Step one: separate out trials by common threads
%   1) target direction
%   2) epoch (BL, AD, WO)
%   3) In the future, task, perturbation type, etc
% Step two: plot trajectories for all trials
%   1) Each dimension against time, color coded by condition
%   2) Any 3 dimensions in trajectory plot, color coded by condition
%
% PARAM_STRUCT OPTIONS:
% Plotting parameters
%   gpfa_params          : parameters of GPFA fits (REQUIRED)
%                            bin_width : (scalar) bin width (in msec) of GPFA model
%                            xdim      : (scalar) number of assumed latent dimensions
%   gpfa_array           : (string) which array to use ('M1','PMd', or 'Both') (REQUIRED)
%   max_trials           : (vector) maximum number of trials to plot of a given condition (Default to all of a given condition)
%   signals              : (cell array) fieldnames of continuous signals to plot (Default to {} for none)
%                            pass in names of fields for continuous signals, e.g. 'vel' or 'acc'
%   trial_conditions     : (cell array) list of trial indices for different conditions. Will plot each element as a different color on the subplots
%                            Can also be vector of trial indices to plot.
%   plot_direction_range : (2 element vector) which range of target directions to plot. Each will get its own columnn. (Defaults to [-pi,pi])
%   plot_dims            : (vector) dimensions to plot against time (Default to first 3 dimensions)
%   align_idx            : (cell array) which idx_EVENT field of trial_data to use for alignment (Defaults to empty, which means start at the beginning of the available data and plot all)
%                            Format is {'idx_EVENT_START',time in msec before; 'idx_EVENT_STOP', time in msec after}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_gpfa(trial_data,params_struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(params_struct,'gpfa_params'), gpfa_params = params_struct.gpfa_params; else gpfa_params = []; end
if isfield(params_struct,'gpfa_array'), gpfa_array = params_struct.gpfa_array; else error('No Array Specified'); end
if isfield(params_struct,'max_trials'), max_trials = params_struct.max_trials; else max_trials = length(trial_data); end
if isfield(params_struct,'signals'), plot_signals = params_struct.signals; else plot_signals = {}; end
if isfield(params_struct,'trial_conditions'), trial_conditions = params_struct.trial_conditions; else trial_conditions = {}; end
if isfield(params_struct,'plot_direction_range'), plot_direction_range = params_struct.plot_direction_range; else plot_direction_range = [-pi,pi]; end
if isfield(params_struct,'plot_dims'), plot_dims = params_struct.plot_dims; else plot_dims = 1:8; end
if isfield(params_struct,'align_idx'), align_idx = params_struct.align_idx; else align_idx = {}; end
if isfield(params_struct,'plot_3d'), plot_3d = params_struct.plot_3d; else plot_3d = false; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   These are parameters that are probably the same all the time
% so it's probably not worth making it an input parameter, but they're here
data_bin_size = 10; %bin size of data in msec
pos_offset = [3-1.7, -33+2.2]; % offset to zero position
event_db = {'idx_target_on','tgt'; ... % list of possible field names for events and a shorthand name
    'idx_go_cue','go'; ...         % add any new events here
    'idx_movement_on','mv'; ...
    'idx_peak_speed','pk'; ...
    'idx_reward','rwd'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   These are a lot of parameters for plotting
% Presumably we won't change these but just in case they are easy to find
font_size          = 12;       % default font size
line_width         = 0.5;        % standard line width
dot_width          = 2;        % standard width of event marker symbols
event_symbol       = 'o';      % standard symbol for event markers
kin_rows           = 3;        % how many rows for kinematic plots
traj_rows          = 4;        % how many rows for time-varying trajectory plots
trial_event_colors = [0    0.4470    0.7410; ... % using default Matlab r2014b color order for trial events
    0.8500    0.3250    0.0980; ...
    0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];
condition_colors  = [0 0 0; ... % colors for different plotting conditions
    1 0 0; ...
    0 0 1; ...
    0 1 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find what directions are within plot_direction_range
all_directions = unique([trial_data.target_direction]);
plot_directions = all_directions( all_directions >= plot_direction_range(1) & all_directions <= plot_direction_range(2) );

if isempty(trial_conditions)
    % default to all trials as one condition
    trial_conditions = {1:length(trial_data)};
elseif ~iscell(trial_conditions)
    % it's probably a single vector of trial indices
    trial_conditions = {trial_conditions};
end

% check for 3d plot dimensions
if plot_3d && length(plot_dims) ~= 3
    error('For 3D plots, please specify only 3 dimensions.');
end

% split out trials based on parameters of trial_conditions and plot_directions
trial_blocks = cell(length(trial_conditions),length(plot_directions));
for iCond = 1:length(trial_conditions)
    for iDir = 1:length(plot_directions)
        trial_blocks{iCond,iDir} = intersect(trial_conditions{iCond}, find(get_trial_data_indices(trial_data,'target_direction',plot_directions(iDir))));
    end
end

% allow for a variable number of events named 'idx_EVENT'
fn = fieldnames(trial_data);
events = fn(cellfun(@(x) ~isempty(regexp(x,'idx_','ONCE')),fn));
clear fn;

% check to see if specified align_idx exists
if ~isempty(align_idx) && ( ~ismember(align_idx{1,1},events) || ~ismember(align_idx{2,1},events) )
    error('Requested alignment index does not exist in trial_data.');
end

if plot_3d
    % find how many rows are needed, this allows for proportional sizing
    num_rows = length(plot_signals)*kin_rows + traj_rows;
    % this is the number of subplots, for use on axis scaling
    num_subplot_rows = length(plot_signals) + 1;
else
    % find how many rows are needed, this allows for proportional sizing
    num_rows = length(plot_signals)*kin_rows + length(plot_dims)*traj_rows;
    % this is the number of subplots, for use on axis scaling
    num_subplot_rows = length(plot_signals) + length(plot_dims);
end

num_cols = length(plot_directions);
% use this to partition the subplot space
subplot_grid = repmat((0:num_rows-1)'*num_cols,1,num_cols) + repmat(1:num_cols,num_rows,1);

% Make new figure
figure('units','normalized','outerposition',[0.1 0 .9 1]);

all_min_x = zeros(1,num_subplot_rows);
all_max_x = zeros(1,num_subplot_rows);
all_min_y = zeros(1,num_subplot_rows);
all_max_y = zeros(1,num_subplot_rows);

% check that data exists and params provided
if isfield(trial_data,[gpfa_array '_gpfa'])
    if ~isempty(gpfa_params)
        gpfa_bin_w = gpfa_params.([gpfa_array '_gpfa']).bin_width;
        gpfa_x_dim = gpfa_params.([gpfa_array '_gpfa']).xdim;
    else
        error('GPFA parameter struct input not provided. See documentation.');
    end
else
    error('GPFA data not present in trial_data.');
end

% loop along conditions
for iCond = 1:length(trial_conditions)
    % loop along desired target directions
    for iDir = 1:length(plot_directions)
        row_count = 1;
        
        temp_trial_data = trial_data(trial_blocks{iCond,iDir});
        num_trials_to_plot = min([length(trial_blocks{iCond,iDir}),max_trials]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop along the specified continuous signals
        % keep count of how many rows have been used
        row_tally = 0;
        for iSignal = 1:length(plot_signals)
            
            % Plot kinematics
            ax(row_count,iDir) = subplot(num_rows,num_cols, ...
                reshape(subplot_grid(row_tally+1:row_tally+kin_rows,iDir)',1,kin_rows ));
            
            if strcmpi(plot_signals{iSignal}(1:3),'pos')
                offset = pos_offset;
            else
                offset = [0 0];
            end
            
            if strcmpi(plot_signals{iSignal}(end-1),'_')
                current_signal = plot_signals{iSignal}(1:end-2);
                % it is likely specified '_x' or '_y'
                switch lower(plot_signals{iSignal}(end-1:end))
                    case '_x'
                        dim_idx = 1;
                    case '_y'
                        dim_idx = 2;
                end
            else
                current_signal = plot_signals{iSignal};
                % plot both
                dim_idx = [1 2];
            end
            
            hold all;
            [min_x,max_x,min_y,max_y] = deal(0);
            for tr_idx = 1:num_trials_to_plot % tr_idx is a dummy variable; useful if you're skipping trials
                if ~isempty(temp_trial_data(tr_idx).(align_idx{1,1})) && ~isempty(temp_trial_data(tr_idx).(align_idx{2,1})) % sometimes an index, especially peak speed, might be missing LOOK INTO
                    % get time indices to plot
                    if isempty(align_idx)
                        t_idx = 1:size(temp_trial_data(tr_idx).(current_signal),1);
                        t_vec = ( 1:size(temp_trial_data(tr_idx).(current_signal),1) ) - 1;
                    else
                        t_start_align = temp_trial_data(tr_idx).(align_idx{1,1});
                        t_stop_align = temp_trial_data(tr_idx).(align_idx{2,1});
                        t_idx = t_start_align - floor(align_idx{1,2}(1)/data_bin_size):t_stop_align + ceil(align_idx{2,2}/data_bin_size);
                        t_vec = ( 1:length(t_idx)) - round(align_idx{1,2}/data_bin_size);
                    end
                    
                    y_vec = temp_trial_data(tr_idx).(current_signal)(t_idx,dim_idx(1)) - offset(dim_idx(1));
                    plot(t_vec,y_vec,'-','LineWidth',line_width,'Color',condition_colors(iCond,:)) % x
                    min_y = min([min_y,min(y_vec)]);
                    max_y = max([max_y,max(y_vec)]);
                    
                    if length(dim_idx) > 1
                        y_vec = temp_trial_data(tr_idx).(current_signal)(t_idx,dim_idx(2)) - offset(dim_idx(2));
                        plot(t_vec,y_vec,'--','LineWidth',line_width,'Color',condition_colors(iCond,:)) % y
                        min_y = min([min_y,min(y_vec)]);
                        max_y = max([max_y,max(y_vec)]);
                    end
                    
                    % keep a tally of max/min values for purposes of setting axes later
                    min_x = min([min_x,t_vec(1)]);
                    max_x = max([max_x,t_vec(end)]);
                end
            end
            ylabel(plot_signals{iSignal},'FontSize',font_size);
            set(gca,'Box','off','TickDir','out','XTick',[],'FontSize',font_size);
            
            % now, add event stamps to each trial
            for tr_idx = 1:num_trials_to_plot
                if ~isempty(temp_trial_data(tr_idx).(align_idx{1,1})) && ~isempty(temp_trial_data(tr_idx).(align_idx{2,1}))  % sometimes an index, especially peak speed, might be missing LOOK INTO
                    % get time indices to plot
                    if isempty(align_idx)
                        t_start_align = 1;
                    else
                        t_start_align = temp_trial_data(tr_idx).(align_idx{1,1});
                    end
                    
                    for iEvent = 1:length(events)
                        plot(temp_trial_data(tr_idx).(events{iEvent})-t_start_align, ...
                            temp_trial_data(tr_idx).(current_signal)(temp_trial_data(tr_idx).(events{iEvent}),dim_idx(1)) - offset(dim_idx(1)), ...
                            event_symbol,'LineWidth',dot_width,'Color',trial_event_colors(iEvent,:));
                        if length(dim_idx) > 1
                            plot(temp_trial_data(tr_idx).(events{iEvent})-t_start_align, ...
                                temp_trial_data(tr_idx).(current_signal)(temp_trial_data(tr_idx).(events{iEvent}),dim_idx(1)) - offset(dim_idx(1)), ...
                                event_symbol,'LineWidth',dot_width,'Color',trial_event_colors(iEvent,:));
                        end
                    end
                end
            end
            
            if row_tally == 0
                title(num2str(plot_directions(iDir)),'FontSize',font_size);
            end
            
            % keep track of the max/min y scaling by subplot row
            all_min_x(row_count) = min([all_min_x(row_count),min_x]);
            all_max_x(row_count) = max([all_max_x(row_count),max_x]);
            all_min_y(row_count) = min([all_min_y(row_count),min_y]);
            all_max_y(row_count) = max([all_max_y(row_count),max_y]);
            
            row_tally = row_tally + kin_rows;
            row_count = row_count + 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot GPFA dimension traces over time
        if plot_3d
            ax(row_count,iDir) = subplot(num_rows,num_cols, ...
                reshape(subplot_grid(row_tally+1:row_tally+traj_rows,iDir)',1,traj_rows ));
            hold all;
            for tr_idx = 1:num_trials_to_plot % tr_idx is a dummy variable; useful if you're skipping trials
                if ~isempty(temp_trial_data(tr_idx).(align_idx{1,1})) && ~isempty(temp_trial_data(tr_idx).(align_idx{2,1}))  % sometimes an index, especially peak speed, might be missing LOOK INTO
                    % get time indices to plot
                    if isempty(align_idx)
                        t_idx = 1:size(temp_trial_data(tr_idx).([gpfa_array '_gpfa']),2);
                    else
                        t_start_align = floor(temp_trial_data(tr_idx).(align_idx{1,1})*(data_bin_size/gpfa_bin_w));
                        t_stop_align = floor(temp_trial_data(tr_idx).(align_idx{2,1})*(data_bin_size/gpfa_bin_w));
                        t_idx = t_start_align - floor(align_idx{1,2}(1)/gpfa_bin_w):t_stop_align + ceil(align_idx{2,2}/gpfa_bin_w);
                    end
                    
                    y_vec = temp_trial_data(tr_idx).([gpfa_array '_gpfa'])(plot_dims,t_idx);
                    plot3(y_vec(1,:),y_vec(2,:),y_vec(3,:),'-','LineWidth',line_width,'Color',condition_colors(iCond,:));
                    
                    % plot event markers
                    for iEvent = 1:length(events)
                        temp_idx = round(temp_trial_data(tr_idx).(events{iEvent})*data_bin_size/gpfa_bin_w);
                        if temp_idx <= t_idx(end) && temp_idx >= t_idx(1)
                            plot3(temp_trial_data(tr_idx).([gpfa_array '_gpfa'])(plot_dims(1), temp_idx ), ...
                                temp_trial_data(tr_idx).([gpfa_array '_gpfa'])(plot_dims(2), temp_idx ), ...
                                temp_trial_data(tr_idx).([gpfa_array '_gpfa'])(plot_dims(3), temp_idx ), ...
                                event_symbol,'LineWidth',dot_width,'Color',trial_event_colors(iEvent,:));
                        end
                    end
                end
            end
            set(gca,'Box','off','TickDir','out','YTick',[],'XTickLabels',[],'FontSize',font_size);
            
            if row_tally == 0
                title(num2str(plot_directions(iDir)),'FontSize',font_size);
            end
            
            row_tally = row_tally + traj_rows;
            row_count = row_count + 1;
        else % do dimensions against time
            for iDim = 1:length(plot_dims)
                if plot_dims(iDim) <= gpfa_x_dim
                    ax(row_count,iDir) = subplot(num_rows,num_cols, ...
                        reshape(subplot_grid(row_tally+1:row_tally+traj_rows,iDir)',1,traj_rows ));
                    hold all;
                    [min_x,max_x,min_y,max_y] = deal(0);
                    for tr_idx = 1:num_trials_to_plot % tr_idx is a dummy variable; useful if you're skipping trials
                        if ~isempty(temp_trial_data(tr_idx).(align_idx{1,1})) && ~isempty(temp_trial_data(tr_idx).(align_idx{2,1}))  % sometimes an index, especially peak speed, might be missing LOOK INTO
                            % get time indices to plot
                            if isempty(align_idx)
                                t_idx = 1:size(temp_trial_data(tr_idx).([gpfa_array '_gpfa']),2);
                                t_vec = ( 1:size(temp_trial_data(tr_idx).([gpfa_array '_gpfa']),2) ) - 1;
                            else
                                t_start_align = floor(temp_trial_data(tr_idx).(align_idx{1,1})*(data_bin_size/gpfa_bin_w));
                                t_stop_align = floor(temp_trial_data(tr_idx).(align_idx{2,1})*(data_bin_size/gpfa_bin_w));
                                t_idx = t_start_align - floor(align_idx{1,2}(1)/gpfa_bin_w):t_stop_align + ceil(align_idx{2,2}/gpfa_bin_w);
                                t_vec = ( 1:length(t_idx)) - floor(align_idx{1,2}/gpfa_bin_w);
                            end
                            
                            y_vec = temp_trial_data(tr_idx).([gpfa_array '_gpfa'])(plot_dims(iDim),t_idx);
                            plot(t_vec,y_vec,'-','LineWidth',line_width,'Color',condition_colors(iCond,:));
                            
                            % keep a tally of max/min values for purposes of setting axes later
                            min_x = min([min_x,t_vec(1)]);
                            max_x = max([max_x,t_vec(end)]);
                            min_y = min([min_y,min(y_vec)]);
                            max_y = max([max_y,max(y_vec)]);
                        end
                    end
                    
                    % plot event markers
                    for tr_idx = 1:num_trials_to_plot
                        if ~isempty(temp_trial_data(tr_idx).(align_idx{1,1})) && ~isempty(temp_trial_data(tr_idx).(align_idx{2,1}))  % sometimes an index, especially peak speed, might be missing LOOK INTO
                            % get time indices to plot
                            if isempty(align_idx)
                                t_start_align = 1;
                            else
                                t_start_align = floor(temp_trial_data(tr_idx).(align_idx{1,1})*(data_bin_size/gpfa_bin_w));
                            end
                            
                            for iEvent = 1:length(events)
                                plot(floor(temp_trial_data(tr_idx).(events{iEvent})*data_bin_size/gpfa_bin_w) - t_start_align, ...
                                    temp_trial_data(tr_idx).([gpfa_array '_gpfa'])(plot_dims(iDim), round(temp_trial_data(tr_idx).(events{iEvent})*data_bin_size/gpfa_bin_w) ), ...
                                    event_symbol,'LineWidth',dot_width,'Color',trial_event_colors(iEvent,:));
                            end
                        end
                    end
                    set(gca,'Box','off','TickDir','out','YTick',[],'XTickLabels',[],'FontSize',font_size);
                    ylabel([gpfa_array ' ' num2str(plot_dims(iDim))],'FontSize',font_size)
                    
                    if row_tally == 0
                        title(num2str(plot_directions(iDir)),'FontSize',font_size);
                    end
                    
                    % again, get the min/max for these rows
                    all_min_x(row_count) = min([all_min_x(row_count),min_x]);
                    all_max_x(row_count) = max([all_max_x(row_count),max_x]);
                    all_min_y(row_count) = min([all_min_y(row_count),min_y]);
                    all_max_y(row_count) = max([all_max_y(row_count),max_y]);
                    
                    row_tally = row_tally + traj_rows;
                    row_count = row_count + 1;
                else
                    warning(['Requested dimension (' num2str(plot_dims(iDim)) ') is larger than available latent dimensions (xDim = ' num2str(gpfa_x_dim) '. Skipping.']);
                end
            end
        end
        xlabel('Time (bins)','FontSize',font_size);
    end
end


if plot_3d
    % link camera rotation for the 3D row
    row_count = row_count - 1;
    Link = linkprop(ax(end,:), {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
    setappdata(gcf, 'StoreTheLink', Link);
end
% go through and set y axes for all plots using the min/max
for iDir = 1:length(plot_directions)
    for iRow = 1:row_count-1
        set(ax(iRow,iDir),'YLim',[all_min_y(iRow), all_max_y(iRow)],'XLim',[all_min_x(iRow), all_max_x(iRow)]);
    end
end

end
