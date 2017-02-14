%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ ] = vis_data( trial_data, params )
%
%   Function to visualize data. Plots in two-column format, where one column
% is a 2-D position plot (and, optionally, a 3-D GPFA trajectory plot), and
% the second column is a stack of time-varying signals (e.g. continuous
% data, spike rasters, GPFA dimensions, etc).
%
% NOTE: needs some tweaking to fix the positions etc
%
% INPUTS:
% Trial_data is a struct array where each element is a trial.
%   target_direction : the angular direction of the target on that trial
%   idx_EVENT        : bin index of trial events. There can be many of these.
%                       Common ones include: target_on, go_cue, movement_on, peak_speed, reward
%   continuous-data  : any number of binned continuous signals (e.g. 'pos','vel','acc','force')
%   ARRAY_spikes     : contains an array [# neurons, # time bins]
%                       Each element is a count of binned spikes. ARRAY is currently 'M1' and/or 'PMd'
%   ARRAY_gpfa       : (if needed) contains an array [# dimensions, # gpfa time bins]
%
% PARAM_STRUCT OPTIONS:
% Plotting parameters
%   trials      : (vector) trial indices to plot. MUST BE SPECIFIED.
%   signals     : (cell array) fieldnames of continuous signals to plot (Default to {'vel,'acc'})
%
% GPFA-specific Parameters
%   plot_gpfa   : (bool) whether to add GPFA data to figure (Default to false, requires GPFA field in trial_data)
%                    NOTE: adds gpfa_dims dimensions to time plots, and adds first 3 dimensions as trajectory plot
%                       To plot only trajectory, pass in empty gpfa_dims parameter
%   gpfa_dims   : (vector) list of GPFA factors to plot
%   gpfa_array  : (string) name of array to plot for GPFA ('M1','PMd', or 'Both')
%   gpfa_params : struct of parameters from run_gpfa (Default to empty)
%
% NOTE: There are a lot more parameters hard-coded at the top of the function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = vis_data( trial_data, params )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
if isfield(params,'trials'), trials_to_plot = params.trials; else, error('No trials specified.'); end
plot_signals      =   {'vel'};
plot_pca          =   false;
pca_dims          =   1:3;
pca_array         =   'M1';
pca_params        =   [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   These are parameters that are less likely to change but can still be
%   overwritten as an input parameter (not documented in help call though)
pos_offset        =   [0, -30]; % offset to zero position
target_size       =   2; % target size in cm
target_distance   =   8; % distance of outer targets in cm
event_db          =   { ...
    'idx_trial_start', 'strt'; ...
    'idx_target_on',   'tgt'; ... % list of possible field names for events and a shorthand name
    'idx_go_cue',      'go'; ...         % add any new events here
    'idx_movement_on', 'mv'; ...
    'idx_peak_speed',  'pk'; ...
    'idx_reward',      'rw'; ...
    'idx_trial_end',   'end'};
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
trial_event_colors =   [0    0.4470    0.7410; ... % using default Matlab r2014b color order for trial events
    0.8500    0.3250    0.0980; ...
    0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin_size     =   trial_data(1).bin_size; %bin size of data in s

% check for foolish inputs
if max(trials_to_plot) > length(trial_data)
    error('Requested too many trials.');
end


num_trials_to_plot = length(trials_to_plot);

% If GPFA plotting is desired, checks that data exists and params provided
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

% add 2D or 3D GPFA trajectory plot below position plot
if plot_pca
    pos_rows_max = floor(num_rows/2);
else
    pos_rows_max = num_rows;
end


% Loop through trials and plot
for tr_idx = 1:num_trials_to_plot % tr_idx is a dummy variable; useful if you're skipping trials
    tr_num = trials_to_plot(tr_idx); % Use tr_num from here down
    
    % Make new figure
    figure('units','normalized','outerposition',[0.1 0 .85 1]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot reach path
    subplot(num_rows,num_cols, ...
        reshape(subplot_grid(1:pos_rows_max,pos_start+1:pos_start+pos_cols)',1,pos_cols*pos_rows_max));
    % above is some logic to partition the available subplot space according to the parameters specified at the top
    hold all;
    plot( ...
        trial_data(tr_num).pos(:,1) - pos_offset(1), ... % x pos
        trial_data(tr_num).pos(:,2) - pos_offset(2), ... % y pos
        'LineWidth',line_width,'Color','k');
    
    % plot center target
    rectangle('Position',[-target_size/2 -target_size/2 target_size target_size],'LineWidth',line_width);
    % plot outer target based on target direction
    x_center = target_distance*cos(trial_data(tr_num).target_direction);
    y_center = target_distance*sin(trial_data(tr_num).target_direction);
    rectangle('Position',[x_center-target_size/2, y_center-target_size/2, target_size, target_size],'LineWidth',line_width);
    axis([pos_range pos_range]);
    axis square;
    set(gca,'XTick',[],'YTick',[], ...
        'YAxisLocation',lower(pos_location),'Box','on','TickDir','out','FontSize',font_size);
    xlabel('X Pos (cm)','FontSize',font_size);
    ylabel('Y Pos (cm)','FontSize',font_size);
    title(['Trial ' num2str(tr_num)],'FontSize',font_size);
    
    % plot dots representing trial events
    for iEvent = 1:length(events)
        plot(trial_data(tr_num).pos(trial_data(tr_num).(events{iEvent}),1) - pos_offset(1), ...
            trial_data(tr_num).pos(trial_data(tr_num).(events{iEvent}),2) - pos_offset(2), ...
            event_symbol,'LineWidth',dot_width,'Color',trial_event_colors(iEvent,:));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot GPFA trajectory path
    if plot_pca
        subplot(num_rows,num_cols, ...
            reshape(subplot_grid(pos_rows_max+1:num_rows,pos_start+1:pos_start+pos_cols)',1,pos_cols*(num_rows-pos_rows_max)));
        hold all;
        plot3(trial_data(tr_num).([pca_array '_pca'])(:,1), ...
            trial_data(tr_num).([pca_array '_pca'])(:,2), ...
            trial_data(tr_num).([pca_array '_pca'])(:,3), ...
            'linewidth',line_width,'color','k');
        
        % now plot event markers
        for iEvent = 1:length(events)
            % scale bin index to fit gpfa bins
            idx = trial_data(tr_num).(events{iEvent});
            
            plot3(trial_data(tr_num).([pca_array '_pca'])(idx,1), ...
                trial_data(tr_num).([pca_array '_pca'])(idx,2), ...
                trial_data(tr_num).([pca_array '_pca'])(idx,3), ...
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
        plot(trial_data(tr_num).(plot_signals{iSignal})(:,1) - offset(1),'r','LineWidth',line_width) % x
        plot(trial_data(tr_num).(plot_signals{iSignal})(:,2) - offset(2),'b','LineWidth',line_width) % y
        xlim([1 size(trial_data(tr_num).pos,1)]);
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
        idx = trial_data(tr_num).(events{iEvent});
        plot([idx idx],[0 1],'LineWidth',line_width,'Color',trial_event_colors(iEvent,:));
        % get the shorthand name from the event_db
        text(idx,.5,event_db{strcmpi(event_db(:,1),events{iEvent}),2},'FontSize',font_size);
    end
    
    axis([1 size(trial_data(tr_num).pos,1) 0 1]);
    set(gca,'Box','off','TickDir','out','XTick',[],'YTick',[],'FontSize',font_size);
    row_tally = row_tally + event_rows;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot PCA dimension traces over time
    if plot_pca
        for iDim = 1:length(pca_dims)
            subplot(num_rows,num_cols, ...
                reshape(subplot_grid(row_tally+1:row_tally+traj_rows,time_start+1:time_start+time_cols)',1,(num_cols-pos_cols)*traj_rows ));
            plot(trial_data(tr_num).([pca_array '_pca'])(:,pca_dims(iDim),:),'k','LineWidth',line_width);
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
        imagesc(trial_data(tr_num).(arrays{iArray})');
        axis([1 size(trial_data(tr_num).pos,1) 0 size(trial_data(tr_num).([arrays{iArray}]),2)]);
        set(gca,'Box','off','TickDir','out','FontSize',font_size);
        ylabel(arrays{iArray}(1:end-7),'FontSize',font_size);
        row_tally = row_tally + spike_rows;
    end
    
    xlabel('Time (bins)','FontSize',font_size);
end

end