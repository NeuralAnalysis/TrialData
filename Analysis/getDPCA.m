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
%     .signals        : (cell) which signals. Two options:
%                           1) {'NAME1','NAME2',etc}
%                           2) {'NAME1',idx; 'NAME2',idx; etc}
%                                   Here idx is which columns to use
%                                   Note: can use 'all' as idx for all
%     .num_dims       : how many dPCA dimensions
%     .do_plot        : flag to make dPCA plot
%
% OUTPUTS:
%
%
% TO DO:
%   - add dPCA projections to trial table
%       Would be getting geometry from trial averaging and projecting
%       single trials (add a guide to say what each dPC is)
%   - support different numbers of decisions for different conditions?
%   - support an arbitrary number of condition structs, rather than one
%   plus the decision_var?
%   - figure out naming parameters and plotting parameters, make them
%   writeable
%
% EXAMPLE:
%   Do dPCA over a few blocks of learning, with target direction as the
%   other factor
% blocks{1} = getTDidx(td,'epoch','BL');
% blocks{2} = getTDidx(td,'epoch','AD','range',[0 0.5]);
% blocks{3} = getTDidx(td,'epoch','AD','range',[0.5 1]);
% getDPCA(td,'target_direction',blocks,struct('signals',{signals}));
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
% Some undocumented extra functions
dpca_plot_fcn  =  @dpca_plot_default;
assignParams(who,params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
if length(conditions) > 2, error('Only two conditions supported at this time.'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THESE SHOULD ALL BECOME INPUT-BASED IN SOME INTUITIVE WAY
% ------------------------------------------------------------------------
% 2. Define parameters
% 1 - condition
% 2 - decision_var
% 3 - time
% [1 3] - condition/time interaction
% [2 3] - decision_var/time interaction
% [1 2] - condition/decision_var interaction
% [1 2 3] - rest
combined_params = { {1,[1 3]}, {2,[2,3]}, {3}, {[1 2],[1 2 3]} };
marg_names      = {'target','learning','time','target/learning interaction'};
marg_colors     = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256; % blue, red, grey, purple
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the size of the massive input matrix
%   N is the number of signals (e.g. neurons)
%   S is the number of conditions --> tasks in our case
%   D is the number of decisions --> targets in our case
%   T is the number of time points --each trial should have the same duration in time !!!
N = sum(cellfun(@(x) length(x),signals(:,2)));
C = cellfun(@length,conditions);
T = size(trial_data(1).(signals{1,1}),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% max number of repetitions
max_trial_num = max(loops_all_the_way_down(0,conditions{:}));

firing_rates        = nan([N,C,T,max_trial_num]);
% Get firing rates
for s = 1:C(1)
    for d = 1:C(2)
        temp = [];
        for n = 1:size(signals,1)
            temp = cat(2,temp,cat(3,trial_data(intersect(conditions{1}{s},conditions{2}{d})).(signals{n,1})));
        end
        for n = 1:size(temp,2)
            firing_rates(n,s,d,:,1:size(temp,3)) = squeeze(temp(:,n,:));
        end
    end
end
% firing_rates_average: N x S x D x T -- these are PSTHs
%   average over the last dimension

firing_rates_avg    = nanmean(firing_rates, 5);

% ------------------------------------------------------------------------
% 2. Do dPCA without regularization

[W, V, which_marg]  = dpca( firing_rates_avg, num_dims, 'combinedParams', combined_params );

expl_var            = dpca_explainedVariance(firing_rates_avg, W, V, ...
    'combinedParams', combined_params);

time_events     = 1;
time            = (1:T)*trial_data(1).bin_size;
if do_plot
    dpca_plot(firing_rates_avg, W, V, dpca_plot_fcn, ...
        'explainedVar', expl_var, ...
        'marginalizationNames', marg_names, ...
        'marginalizationColours', marg_colors, ...
        'whichMarg', which_marg,                 ...
        'time', time,                        ...
        'timeEvents', time_events,               ...
        'timeMarginalization', 3, ...
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
% I'm pretty proud of this one. Will find the largest number of common
% trials for any arbitrary number of condition inputs, using recursive
% loops
function val = loops_all_the_way_down(val,varargin)
% get all of the trials that are common between the conditions
if all(cellfun(@isnumeric,varargin))
    % we've reached the bottom of the hole. Find the common trials
    for i = 1:length(varargin)-1
        val = [val,length(intersect(varargin{i:i+1}))];
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
function [val,valIdx] = loops_all_the_way_down2(depthCount,valIdx,val,trial_data,signals,max_trial_num,varargin)
if all(cellfun(@isnumeric,varargin))
    valIdx = [valIdx, depthCount];
    % we've reached the bottom of the hole. Start computing
    % NOTE currently hard coded assuming that each condition cell will be a
    % 3-D array (one D for time, one D for trials, and one D for the
    % condition variable)
    fr = [];
    % concatenate together all data for all of the requested signals
    all_trials = 1:length(trial_data);
    for i = 1:length(varargin)-1
        all_trials = intersect(all_trials,varargin{i});
    end
    for n = 1:size(signals,1)
        fr = cat(2,fr,cat(3,trial_data(all_trials).(signals{n,1})));
    end
    fr = permute(fr,[2,1,3]);
    temp = NaN(sum(cellfun(@length,signals(:,2))),size(trial_data(1).(signals{1,1}),1),max_trial_num);
    temp(:,:,1:size(fr,3)) = fr;
    val = [val {temp}];
else
    % still need to loop further. Find the next thing that is still a
    % cell and start looping along it
    cell_idx = find(cellfun(@iscell,varargin),1,'first');
    iter_data = varargin{cell_idx};
    for i = 1:length(iter_data)
        temp = varargin;
        temp{cell_idx} = iter_data{i};
        [val,valIdx] = loops_all_the_way_down2([depthCount;i],valIdx,val,trial_data,signals,max_trial_num,temp{:});
    end
end

% check to see if this is the end of the infinite loop
if all(cellfun(@iscell,varargin)) && length(val) == prod(cellfun(@(x) length(cellfun(@numel,x)),varargin))
    % now process everything and concatenate together
    num_conds = size(valIdx,1);
    idx = cell(1,num_conds);
    for i = 1:num_conds
        idx_vals{i} = unique(valIdx(i,:));
    end
    val = loop_me_twice(length(size(val{1})),[],val,valIdx,0,idx_vals{:});
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A recursive sub-function for my recursive fuction
function valCat = loop_me_twice(dim,valCat,val,valIdx,depthCount,varargin)
depthCount=depthCount+1;
if length(varargin) == 1
        temp = cat(dim,val{:});%{idx(x,:)==varargin{1}(i)});
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
