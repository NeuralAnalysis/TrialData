%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data, dPCA_info] = getDPCA(trial_data,varargin)
%
%   Finds dPCA for a dataset according to arbitrary conditions.
%
%   [trial_data, dPCA_info] = getDPCA(trial_data, ...conditions... params)
%       Input any number of conditions. The last input should always be
% params struct. Conditions are strings (for fields, e.g. 'target_direction')
% or cell arrays of trial indices.
%
%
% INPUTS:
%   params : parameter struct
%     .signals         : (cell) which signals. Two options:
%                           1) {'NAME1','NAME2',etc}
%                           2) {'NAME1',idx; 'NAME2',idx; etc}
%                                   Here idx is which columns to use
%                                   Note: can use 'all' as idx for all
%     .num_dims        : how many dPCA dimensions
%     .do_plot         : flag to make dPCA plot
%     .combined_params : see dPCA code or example below
%     .marg_names      : name of each marginalization
%     .marg_colors     : color for each marginalization
%
% OUTPUTS:
%   trial_data : data struct with dPCA projections added
%   dPCA_info  : struct with info about dPCA analysis
%
% TO DO:
%   - add dPCA projections to trial table
%       Would be getting geometry from trial averaging and projecting
%       single trials (add a guide to say what each dPC is)
%   - figure out naming parameters and plotting parameters, make them
%        writeable
%   - support different numbers of conditions for different conditions?
%           e.g. if Condition1 is task and Condition2 is target, can have 6
%           targets in one task and 8 targets in another. Currently you
%           need to have the same number of targets for all tasks
%
% Written by Matt Perich. Adapted from Juan's code. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data, dPCA_info] = getDPCA(trial_data,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = varargin{end}; if ~isstruct(params), error('Last input must be params.'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
signals        =  [];
num_dims       =  15;
do_plot        =  true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are some more complicated inputs for the dPCA code
% 1 - condition
% 2 - decision_var
% 3 - time
% [1 3] - condition/time interaction
% [2 3] - decision_var/time interaction
% [1 2] - condition/decision_var interaction
% [1 2 3] - rest
combined_params = { {1}, {2,[1 2]}, {3,[1,3]}, {[2 3],[1 2 3]} };
marg_names      = {'time','target','learning','target/learning interaction'};
marg_colors     = [150 150 150; 23 100 171; 187 20 25; 114 97 171]/256; % blue, red, grey, purple
W = [];
V = [];
which_marg = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra functions
dpca_plot_fcn  =  @dpca_plot_td;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assignParams(who,params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end
signals        = check_signals(trial_data(1),signals);
% make sure all trials are same length in time
if length(unique(cellfun(@(x) size(x,1),{trial_data.(signals{1,1})}))) > 1
    error('Trials are not all the same length.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse out condition inputs, build a cell array of cell arrays
conditions = cell(1,length(varargin)-1);
for i = 1:length(varargin)-1
    if ischar(varargin{i}) % it's a fieldname
        if ~isfield(trial_data,varargin{i}), error('Condition input fieldname not recognized.'); end
        % get unique values for this field name
        if ischar(trial_data(1).(varargin{i}))
            u = unique({trial_data.(varargin{i})});
        else
            u = num2cell(unique([trial_data.(varargin{i})]));
        end
        % now get the indices for each trial of each type
        conditions{i} = cellfun(@(x) getTDidx(trial_data,varargin{i},x),u,'uni',0);
    elseif iscell(varargin{i}) % it's a group of trials
        conditions{i} = varargin{i};
    else
        error('Condition input should be string fieldname or cell array of trial indices.');
    end
end
if length(conditions) > 3, warning('This many conditions takes a while to run...'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the number of time points
T = size(trial_data(1).(signals{1,1}),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% max number of repetitions
max_trial_num = max(loops_all_the_way_down(0,conditions{:}));

% firing_rates: N x T x max_trial_num x cond1 x cond2 x ... etc
firing_rates = loops_for_fr([],[],{},trial_data,signals,max_trial_num,conditions{:});
% The order will always be neurons first, time second, then the other
% conditions in the OPPOSITE order of conditions. So I flip those
firing_rates = permute(firing_rates,[1:2,2+fliplr(1:length(conditions))]);

% ------------------------------------------------------------------------
% 2. Do dPCA without regularization

if isempty(W)
    [W, V, which_marg]  = dpca( firing_rates, num_dims, 'combinedParams', combined_params );
end

expl_var            = dpca_explainedVariance(firing_rates, W, V, ...
    'combinedParams', combined_params);

time_events     = 1;
time            = (1:T)*trial_data(1).bin_size;
if do_plot
    dpca_plot(firing_rates, W, V, dpca_plot_fcn, ...
        'explainedVar', expl_var, ...
        'marginalizationNames', marg_names, ...
        'marginalizationColours', marg_colors, ...
        'whichMarg', which_marg,                 ...
        'time', time,                        ...
        'timeEvents', time_events,               ...
        'timeMarginalization', 1, ...
        'legendSubplot', num_dims);
end


% ------------------------------------------------------------------------
% 2. Do dPCA in each marginalization separately

% dpca_perMarginalization(firing_rates_avg, @dpca_plot_default, ...
%    'combinedParams', combined_params);


% ------------------------------------------------------------------------
% project data onto dPC axes
% [lat_vars, lat_vars_st] = get_lat_vars_dPCA( firing_rates_avg, firing_rates, W, ...
%     'explainedVar', expl_var, ...
%     'marginalizationNames', marg_names, ...
%     'marginalizationColours', marg_colors, ...
%     'whichMarg', which_marg,                 ...
%     'time', time,                        ...
%     'timeEvents', time_events,               ...
%     'timeMarginalization', 3, ...
%     'legendSubplot', num_dims);


% ------------------------------------------------------------------------
% 3. Return dPCA results

dPCA_info.W               = W;
dPCA_info.V               = V;
% dPCA_info.lat_vars_mn     = lat_vars;
% dPCA_info.lat_vars        = lat_vars_st;
dPCA_info.which_marg      = which_marg;
dPCA_info.marg_names      = marg_names;
dPCA_info.expl_var        = expl_var;
dPCA_info.num_comps       = num_dims;
dPCA_info.combined_params = combined_params;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Will find the largest number of common trials for any arbitrary number
% of condition inputs, using recursive loops
function val = loops_all_the_way_down(val,varargin)
    % get all of the trials that are common between the conditions
    if all(cellfun(@isnumeric,varargin))
        % we've reached the bottom of the hole. Find the common trials
        if length(varargin) == 1
            val = [val,length(varargin{1})];
        else
            for i = 1:length(varargin)-1
                val = [val,length(intersect(varargin{i:i+1}))];
            end
        end
    else
        % still need to loop further. Find the next thing that is still a
        % cell and start looping along it
        cell_idx = find(cellfun(@iscell,varargin),1,'first');
        iter_data = varargin{cell_idx};
        for i = 1:length(iter_data)
            temp = varargin;
            temp{cell_idx} = iter_data{i};
            val = loops_all_the_way_down(val,temp{:});
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Another crazy recursive function to split trial_data up by trials
function [val,valIdx] = loops_for_fr(depthCount,valIdx,val,trial_data,signals,max_trial_num,varargin)
    if all(cellfun(@isnumeric,varargin))
        valIdx = [valIdx, depthCount];
        % we've reached the bottom of the hole. Start computing
        % NOTE currently hard coded assuming that each condition cell will be a
        % 2-D array (one D for time and one D for the condition variable)
        fr = [];
        % find all common trials
        all_trials = varargin{1};
        for i = 2:length(varargin)
            all_trials = intersect(all_trials,varargin{i});
        end
        for n = 1:size(signals,1)
            temp_fr = cat(3,trial_data(all_trials).(signals{n,1}));
            temp_fr = temp_fr(:,signals{n,2},:);
            fr = cat(2,fr,temp_fr);
        end
        % if no trial matches, this will crash
        if ~isempty(fr), fr = mean(fr,3); end
        val = [val {fr'}];
    else
        % still need to loop further. Find the next thing that is still a
        % cell and start looping along it
        cell_idx = find(cellfun(@iscell,varargin),1,'first');
        iter_data = varargin{cell_idx};
        for i = 1:length(iter_data)
            temp = varargin;
            temp{cell_idx} = iter_data{i};
            [val,valIdx] = loops_for_fr([depthCount;i],valIdx,val,trial_data,signals,max_trial_num,temp{:});
        end
    end
    
    % check to see if this is the end of the infinite loop
    if all(cellfun(@iscell,varargin)) && length(val) == prod(cellfun(@(x) length(cellfun(@numel,x)),varargin))
        % now process everything and concatenate together
        num_conds = size(valIdx,1);
        idx_vals = cell(1,num_conds);
        for i = 1:num_conds
            idx_vals{i} = unique(valIdx(i,:));
        end
        if length(varargin) == 1 % only one input is easy
            val = cat(3,val{:});
        else
            val = loop_me_twice(length(size(val{1})),[],val,valIdx,0,idx_vals{:});
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A recursive sub-function for my recursive fuction
function valCat = loop_me_twice(dim,valCat,val,valIdx,depthCount,varargin)
    depthCount=depthCount+1;
    if length(varargin) == 1
            temp = cat(dim,val{:});
            valCat = cat(dim+1,temp,valCat);
    else
        for i = 1:length(varargin{1})
            p = valIdx(depthCount,:)==i;
            valCat = loop_me_twice(dim+1,valCat,val(p),valIdx(:,p),depthCount,varargin{2:end});
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
