%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [avg_data, cond_idx] = trialAverage(trial_data, params)
%
%   Returns struct where each entry is an average for all time-varying
% signals across trials for some condition. Requires every trial be the
% same length. This can be easily done with trimTD or stretchSignals.
%
% Note: returns all fields, but for meta parameters, etc, you can average
% across trials with different values. If they all are the same, fills in
% field for that condition with that value. Otherwise, adds cell array that
% contains all unique values, largely for later reference.
%
% INPUTS:
%   trial_data : the struct
%   params     : struct with parameters
%     .conditions : (string or cell array) the field name(s) in trial_data
%                       avg_data will have an entry for each unique combo
%                       Hint: use 'all' to just do all of trial_data
%     .add_std    : (bool) whether to add a field with the standard dev at
%                       each time point. it will be trial_data.SIGNALNAME_std
%
%    NOTE: can pass conditions (string/cell) as second input instead of
%          params struct.
%
% OUTPUTS:
%   avg_data   : struct representing average across trials for each condition
%   cond_idx   : cell array containing indices for each condition
%
% EXAMPLES:
%   e.g. to average over all target directions and task epochs
%       avg_data = trialAverage(trial_data,{'target_direction','epoch'});
%       Note: gives a struct of size #_TARGETS * #_EPOCHS
%
% TO DO:
%   Add a check that unit_guides and lfp_guides are the same for all trials
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [avg_data,cond_idx] = trialAverage(trial_data, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
add_std     =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
conditions  =  {};
verbose = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2 && isstruct(params) % overwrite parameters
    assignParams(who,params);
elseif nargin == 2 % conditions is passed in
    conditions = params;
elseif nargin < 2 % do all
    conditions = {'all'};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
if isempty(conditions), conditions = 'all'; end
if strcmpi(conditions,'all')
    disp('trialAverage: No conditions provided. Averaging all!');
end
if ~iscell(conditions), conditions = {conditions}; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get list of time-varying signals that we will average over
time_vars = getTDfields(trial_data,'time');

fn_time = getTDfields(trial_data,'time');
if length(unique(cellfun(@(x) size(x,1),{trial_data.(fn_time{1})}))) ~= 1
    error('Trials are not uniform length. Do this with time stretching option or trimTD.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(conditions{1},'all')
    all_conds = {'all'};
    num_conds = 1;
else
    % loop along conditions and get unique values for each
    cond_vals = cell(1,length(conditions));
    for iCond = 1:length(conditions)
        % unique doesn't work on numeric cell arrays for some reason
        if ischar(trial_data(1).(conditions{iCond}))
            uc = unique({trial_data.(conditions{iCond})});
        elseif isnumeric(trial_data(1).(conditions{iCond})) || islogical(trial_data(1).(conditions{iCond}))
            uc = num2cell(unique([trial_data.(conditions{iCond})]));
        end
        cond_vals{iCond} = uc;
    end
    % build a list of all possible combinations of values
    temp = cellfun(@(x) 1:length(x),cond_vals,'Uni',0);
    all_conds = combvec(temp{:})';
    num_conds = size(all_conds,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get indices for trials meeting each unique combination
%   and make struct with trial-averaged entries
fn_meta         =  getTDfields(trial_data,'meta');
fn_unit_guides  =  getTDfields(trial_data,'unit_guides');
fn_lfp_guides   =  getTDfields(trial_data,'lfp_guides');

cond_idx = cell(1,num_conds);
avg_data = repmat(struct(),1,num_conds);
for trial = 1:num_conds
    if strcmpi(conditions{1},'all')
        cond_idx{trial} = 1:length(trial_data);
    else
        func_in = cell(1,2*size(all_conds,2));
        for iCond = 1:size(all_conds,2)
            func_in{2*(iCond-1)+1}   = conditions{iCond};
            func_in{2*(iCond-1)+2} = cond_vals{iCond}{all_conds(trial,iCond)};
            avg_data(trial).(conditions{iCond}) = cond_vals{iCond}{all_conds(trial,iCond)};
        end
        cond_idx{trial}=getTDidx(trial_data,func_in);
    end
    
    % populate meta fields
    for f = 1:length(fn_meta)
        if ~isempty(regexp(fn_meta{f},'_names', 'once')) % names are a special kind of meta
            avg_data(trial).(fn_meta{f}) = trial_data(cond_idx{trial}).(fn_meta{f});
        elseif ischar(trial_data(1).(fn_meta{f}))
            u = unique({trial_data(cond_idx{trial}).(fn_meta{f})});
            if length(u) == 1
                avg_data(trial).(fn_meta{f}) = u{1};
            else
                avg_data(trial).(fn_meta{f}) = u;
            end
        elseif iscell(trial_data(1).(fn_meta{f}))
            avg_data(trial).(fn_meta{f}) = unique([trial_data(cond_idx{trial}).(fn_meta{f})]);
        else
            if size(trial_data(1).(fn_meta{f}),2) > 1
                u = unique(cat(1,trial_data(cond_idx{trial}).(fn_meta{f})),'rows');
            else
                u = unique([trial_data(cond_idx{trial}).(fn_meta{f})]);
            end
            avg_data(trial).(fn_meta{f}) = u;
        end
    end
    
    % add unit guides
    for f = 1:length(fn_unit_guides)
        avg_data(trial).(fn_unit_guides{f}) = trial_data(1).(fn_unit_guides{f});
    end
    
    % add lfp guides
    for f = 1:length(fn_lfp_guides)
        avg_data(trial).(fn_lfp_guides{f}) = trial_data(1).(fn_lfp_guides{f});
    end
    
    % now loop along time signals to average
    for v = 1:length(time_vars)
        avg_data(trial).(time_vars{v}) = mean(cat(3,trial_data(cond_idx{trial}).(time_vars{v})),3);
        if add_std
            avg_data(trial).([time_vars{v} '_std']) = std(cat(3,trial_data(cond_idx{trial}).(time_vars{v})),[],3);
        end
    end
    
    % add idx fields
    fn_idx = getTDfields(trial_data,'idx');
    
    for f = 1:length(fn_idx)
        if length(unique([trial_data(cond_idx{trial}).(fn_idx{f})])) == 1
            avg_data(trial).(fn_idx{f}) = trial_data(cond_idx{trial}(1)).(fn_idx{f});
        end
    end
    
end

% restore logical order
avg_data = reorderTDfields(avg_data);


