%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ ] = visData( trial_data, params )
%
%   Function to visualize data. Plots in two-column format, where one column
% is a 2-D position plot (and, optionally, a 3-D PCA trajectory plot), and
% the second column is a stack of time-varying signals (e.g. continuous
% data, spike rasters, PCA dimensions, etc).
%
% NOTE: needs some tweaking to fix the positions etc
%
% INPUTS:
%   Trial_data : struct array where each element is a trial.
%
%   params     : parameter struct
%     .trials      : (vector) trial indices to plot
%                       Note: must be specified unless a single trial is passed in
%     .signals     : (cell array) fieldnames of continuous signals to plot (Default to {'vel,'acc'})
%     .target_direction : the angular direction of the target on that trial
%     .idx_EVENT        : bin index of trial events. There can be many of these.
%                       Common ones include: target_on, go_cue, movement_on, peak_speed, reward
%     .continuous-data  : any number of binned continuous signals (e.g. 'pos','vel','acc','force')
%     .ARRAY_spikes     : contains an array [# neurons, # time bins]
%                       Each element is a count of binned spikes. ARRAY is currently 'M1' and/or 'PMd'
%     .ARRAY_pca       : (if needed) contains an array [# dimensions, # pca time bins]
%
% PCA-specific Parameters (should generalize this functionality more later)
%     .plot_pca   : (bool) whether to add PCA data to figure (Default to false, requires _pca field in trial_data)
%                    NOTE: adds pca_dims dimensions to time plots, and adds first 3 dimensions as trajectory plot
%                       To plot only trajectory, pass in empty pca_dims parameter
%     .pca_dims   : (vector) list of PCA dims to plot
%     .pca_array  : (string) name of array to plot for PCA ('M1','PMd', or 'Both')
% NOTE: There are a lot more parameters hard-coded at the top of the function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = visData( trial_data, params )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
trials_to_plot    =   1;
plot_signals      =   {'vel'};
plot_pca          =   false;
pca_dims          =   1:3;
pca_array         =   '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   These are parameters that are less likely to change but can still be
%   overwritten as an input parameter (not documented in help call though)
pos_offset        =   [0, 0]; % offset to zero position
target_size       =   2; % target size in cm
target_distance   =   8; % distance of outer targets in cm
event_db          =   { ...
    'idx_trial_start', 'strt'; ...
    'idx_startTime',   'strt'; ...
    'idx_target_on',   'tgt'; ... % list of possible field names for events and a shorthand name
    'idx_tgtOnTime',   'tgt'; ...         % add any new events here
    'idx_bumpTime',    'bump'; ...
    'idx_go_cue',      'go'; ...
    'idx_goCueTime',   'go'; ...
    'idx_movement_on', 'mv'; ...
    'idx_peak_speed',  'pk'; ...
    'idx_reward',      'rw'; ...
    'idx_trial_end',   'end'; ...
    'idx_endTime',     'end'};
verbose           =    false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are a lot of parameters for plotting
%   Presumably we won't change these but just in case you can
pos_range          =   [-9,9];   % range for 2-D position plot axes
font_size          =   12;       % default font size
line_width         =   2;        % standard line width
dot_width          =   3;        % standard width for dot markers
event_symbol       =   'o';      % standard symbol for event markers
pos_cols           =   2;        % how many columns for position plot
time_cols          =   3;        % how many columns for time-variable plots
kin_rows           =   3;        % how many rows for kinematic plots
event_rows         =   1;        % how many rows for event markers
traj_rows          =   3;        % how many rows for time-varying trajectory plots
spike_rows         =   4;        % how many rows for spike markers
pos_location       =   'right'; % if position plot is on 'left' or 'right'
trial_event_colors =   parula(size(event_db,1)); % use default matlab colors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1, assignParams(who,params); else params = struct(); end% overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bin_size     =   trial_data(1).bin_size; %bin size of data in s
if ~isfield(params,'trials_to_plot') && length(trial_data) > 1
    error('No trials specified.');
end
% check for foolish inputs
if max(trials_to_plot) > length(trial_data)
    error('Requested too many trials.');
end

if isempty(pca_array) && plot_pca
    disp('WARNING: No PCA array provided. Skipping PCA plots.');
    plot_pca = false;
end

num_trials_to_plot = length(trials_to_plot);

% If PCA plotting is desired, checks that data exists and params provided
if plot_pca && ~isfield(trial_data,[pca_array '_pca'])
    error('PCA data not present in trial_data.');
end

% allow for a variable number of arrays and events
%   assumes each array has a field named ARRAY_spikes
fn = fieldnames(trial_data);
arrays = fn(cellfun(@(x) ~isempty(regexp(x,'_spikes','ONCE')),fn));
num_arrays = length(arrays);
events = fn(cellfun(@(x) ~isempty(regexp(x,'idx_','ONCE')),fn));
clear fn;

% find how many rows are needed
num_rows = length(plot_signals)*kin_rows + event_rows + num_arrays*spike_rows;
if plot_pca, num_rows=num_rows+length(pca_dims)*traj_rows; end

num_cols = pos_cols + time_cols;
% use this to partition the subplot space
subplot_grid = repmat((0:num_rows-1)'*num_cols,1,num_cols) + repmat(1:num_cols,num_rows,1);

% some variables to position the columns
switch lower(pos_location)
    case 'left'
        pos_start = 0;
        time_start = pos_cols;
    case 'right'
        time_start = 0;
        pos_start = time_cols;
end

% add 2D or 3D PCA trajectory plot below position plot
if plot_pca
    pos_rows_max = floor(num_rows/2);
else
    pos_rows_max = num_rows;
end


% Loop through trials and plot
for tr_idx = 1:num_trials_to_plot % tr_idx is a dummy variable; useful if you're skipping trials
    trial = trials_to_plot(tr_idx); % Use tr_num from here down
    
    % check to make sure events aren't empty or nan
    idx = true(1,length(events));
    for iEvent = 1:length(events)
        if isempty(trial_data(trial).(events{iEvent})) || isnan(trial_data(trial).(events{iEvent})) 
            idx(iEvent) = false;
        end
    end
    events = events(idx); clear idx;
    
    % Make new figure
    figure('units','normalized','outerposition',[0.1 0 .85 1]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot reach path
    subplot(num_rows,num_cols, ...
        reshape(subplot_grid(1:pos_rows_max,pos_start+1:pos_start+pos_cols)',1,pos_cols*pos_rows_max));
    % above is some logic to partition the available subplot space according to the parameters specified at the top
    hold all;
    plot( ...
        trial_data(trial).pos(:,1) - pos_offset(1), ... % x pos
        trial_data(trial).pos(:,2) - pos_offset(2), ... % y pos
        'LineWidth',line_width,'Color','k');
    
    % plot center target
    rectangle('Position',[-target_size/2 -target_size/2 target_size target_size],'LineWidth',line_width);
    % plot outer target based on target direction
    x_center = target_distance*cos(trial_data(trial).target_direction);
    y_center = target_distance*sin(trial_data(trial).target_direction);
    rectangle('Position',[x_center-target_size/2, y_center-target_size/2, target_size, target_size],'LineWidth',line_width);
    axis([pos_range pos_range]);
    axis square;
    set(gca,'XTick',[],'YTick',[], ...
        'YAxisLocation',lower(pos_location),'Box','on','TickDir','out','FontSize',font_size);
    xlabel('X Pos (cm)','FontSize',font_size);
    ylabel('Y Pos (cm)','FontSize',font_size);
    title(['Trial ' num2str(trial)],'FontSize',font_size);
    
    % plot dots representing trial events
    for iEvent = 1:length(events)
        plot(trial_data(trial).pos(trial_data(trial).(events{iEvent}),1) - pos_offset(1), ...
            trial_data(trial).pos(trial_data(trial).(events{iEvent}),2) - pos_offset(2), ...
            event_symbol,'LineWidth',dot_width,'Color',trial_event_colors(iEvent,:));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot PCA trajectory path
    if plot_pca
        subplot(num_rows,num_cols, ...
            reshape(subplot_grid(pos_rows_max+1:num_rows,pos_start+1:pos_start+pos_cols)',1,pos_cols*(num_rows-pos_rows_max)));
        hold all;
        plot3(trial_data(trial).([pca_array '_pca'])(:,1), ...
            trial_data(trial).([pca_array '_pca'])(:,2), ...
            trial_data(trial).([pca_array '_pca'])(:,3), ...
            'linewidth',line_width,'color','k');
        
        % now plot event markers
        for iEvent = 1:length(events)
            % get indices
            idx = trial_data(trial).(events{iEvent});
            
            plot3(trial_data(trial).([pca_array '_pca'])(idx,1), ...
                trial_data(trial).([pca_array '_pca'])(idx,2), ...
                trial_data(trial).([pca_array '_pca'])(idx,3), ...
                event_symbol,'linewidth',dot_width,'color',trial_event_colors(iEvent,:));
        end
        
        set(gca,'XTick',[],'YTick',[],'ZTick',[], ...
            'Box','off','TickDir','out','FontSize',font_size);
        xlabel('PC1','FontSize',font_size);
        ylabel('PC2','FontSize',font_size);
        zlabel('PC3','FontSize',font_size);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop along the specified continuous signals
    % keep count of how many rows have been used
    row_tally = 0;
    for iSignal = 1:length(plot_signals)
        % Plot kinematics
        subplot(num_rows,num_cols, ...
            reshape(subplot_grid(row_tally+1:row_tally+kin_rows,time_start+1:time_start+time_cols)',1,(num_cols-pos_cols)*kin_rows ));
        
        if strcmpi(plot_signals{iSignal},'pos')
            offset = pos_offset;
        else
            offset = [0 0];
        end
        
        hold all;
        plot(trial_data(trial).(plot_signals{iSignal})(:,1) - offset(1),'r','LineWidth',line_width) % x
        plot(trial_data(trial).(plot_signals{iSignal})(:,2) - offset(2),'b','LineWidth',line_width) % y
        xlim([1 size(trial_data(trial).pos,1)]);
        ylabel(plot_signals{iSignal},'FontSize',font_size);
        set(gca,'Box','off','TickDir','out','XTick',[],'FontSize',font_size);
        row_tally = row_tally + kin_rows;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot events
    subplot(num_rows,num_cols, ...
        reshape(subplot_grid(row_tally+1:row_tally+event_rows,time_start+1:time_start+time_cols)',1,(num_cols-pos_cols)*event_rows ));
    hold all;
    
    for iEvent = 1:length(events)
        idx = trial_data(trial).(events{iEvent});
        plot([idx idx],[0 1],'LineWidth',line_width,'Color',trial_event_colors(iEvent,:));
        % get the shorthand name from the event_db
        if any(strcmpi(event_db(:,1),events{iEvent}))
            text(idx,.5,event_db{strcmpi(event_db(:,1),events{iEvent}),2},'FontSize',font_size);
        else
            text(idx,.5,strrep(events{iEvent},'idx_',''),'FontSize',font_size);
        end
    end
    
    axis([1 size(trial_data(trial).pos,1) 0 1]);
    set(gca,'Box','off','TickDir','out','XTick',[],'YTick',[],'FontSize',font_size);
    row_tally = row_tally + event_rows;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot PCA dimension traces over time
    if plot_pca
        for iDim = 1:length(pca_dims)
            subplot(num_rows,num_cols, ...
                reshape(subplot_grid(row_tally+1:row_tally+traj_rows,time_start+1:time_start+time_cols)',1,(num_cols-pos_cols)*traj_rows ));
            plot(trial_data(trial).([pca_array '_pca'])(:,pca_dims(iDim),:),'k','LineWidth',line_width);
            axis('tight');
            set(gca,'Box','off','TickDir','out','YTick',[],'XTickLabels',[],'FontSize',font_size);
            ylabel([pca_array ' PC' num2str(pca_dims(iDim))],'FontSize',font_size)
            
            row_tally = row_tally + traj_rows;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot spiking rasters
    for iArray = 1:num_arrays
        set(gca,'XTickLabels',[]);
        subplot(num_rows,num_cols, ...
            reshape(subplot_grid(row_tally+1:row_tally+spike_rows,time_start+1:time_start+time_cols)',1,(num_cols-pos_cols)*spike_rows ));
        imagesc(trial_data(trial).(arrays{iArray})');
        axis([1 size(trial_data(trial).pos,1) 0 size(trial_data(trial).([arrays{iArray}]),2)]);
        set(gca,'Box','off','TickDir','out','FontSize',font_size);
        ylabel(arrays{iArray}(1:end-7),'FontSize',font_size);
        row_tally = row_tally + spike_rows;
    end
    
    xlabel('Time (bins)','FontSize',font_size);
end

end
