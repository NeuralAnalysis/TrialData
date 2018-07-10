%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [avg_data, cond_idx] = trialAverage(trial_data, conditions, params)
%
%   Returns struct where each entry is an average for all time-varying
% signals across trials for some condition. This function can interpolate
% to stretch/shrink all trials to the same number of points, but it does
% not trim data. This can be easily done with trimTD.
%
% Note: returns all fields, but for meta parameters, etc, you can average
% across trials with different values. If they all are the same, fills in
% field for that condition with that value. Otherwise, adds cell array that
% contains all unique values, largely for later reference.
%
% INPUTS:
%   trial_data : the struct
%   conditions : (string or cell array) the field name(s) in trial_data
%                       avg_data will have an entry for each unique combo
%                       Hint: use 'all' to just do all of trial_data
%   params     : struct with parameters
%     .do_stretch : (bool) whether to stretch/shrink trials to uniform length
%                       if false, all trials must have same number of points
%     .num_samp   : (int) how many time points to use for interpolation
%     .add_std    : (bool) whether to add a field with the standard dev at
%                       each time point. it will be trial_data.SIGNALNAME_std
%
% OUTPUTS:
%   avg_data : struct representing average across trials for each condition
%   cond_idx : cell array containing indices for each condition
%
% EXAMPLES:
%   e.g. to average over all target directions and task epochs
%       avg_data = trialAverage(trial_data,{'target_direction','epoch'});
%       Note: gives a struct of size #_TARGETS * #_EPOCHS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [avg_data,cond_idx] = trialAverage(trial_data, conditions, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
do_stretch  =  false;
num_samp    =  1000;
add_std     =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
avg_flag = true; % will add a flag field saying it's trial-averaged
if nargin > 2, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end
if nargin == 1, error('Conditions not provided as input.'); end
if ~iscell(conditions), conditions = {conditions}; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get list of time-varying signals that we will average over
time_vars = getTDfields(trial_data,'time');

if strcmpi(conditions,'all'), error('all not implemented yet...'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time warp each trial to the same number of points, if desired
if do_stretch
    trial_data = stretchSignals(trial_data,struct('num_samp',num_samp));
end

fn_time = getTDfields(trial_data,'time');
if length(unique(cellfun(@(x) size(x,1),{trial_data.(fn_time{1})}))) ~= 1
    error('Trials are not uniform length. Do this with time stretching option or trimTD.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
temp=cellfun(@(x) 1:length(x),cond_vals,'Uni',0);
all_conds = combvec(temp{:})';
num_conds = size(all_conds,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get indices for trials meeting each unique combination
%   and make struct with trial-averaged entries
fn_meta = getTDfields(trial_data,'meta');

cond_idx = cell(1,num_conds);
avg_data = repmat(struct(),1,num_conds);
for i = 1:num_conds
    func_in = cell(1,2*size(all_conds,2));
    for iCond = 1:size(all_conds,2)
        func_in{2*(iCond-1)+1}   = conditions{iCond};
        func_in{2*(iCond-1)+2} = cond_vals{iCond}{all_conds(i,iCond)};
        avg_data(i).(conditions{iCond}) = cond_vals{iCond}{all_conds(i,iCond)};
    end
    cond_idx{i}=getTDidx(trial_data,func_in);
    
    % populate meta fields
    for f = 1:length(fn_meta)
        if ~isempty(regexp(fn_meta{f},'_names', 'once')) % names are a special kind of meta
            avg_data(i).(fn_meta{f}) = trial_data(cond_idx{i}).(fn_meta{f});
        elseif ischar(trial_data(1).(fn_meta{f}))
            u = unique({trial_data(cond_idx{i}).(fn_meta{f})});
            if length(u) == 1
                avg_data(i).(fn_meta{f}) = u{1};
            else
                avg_data(i).(fn_meta{f}) = u;
            end
        elseif iscell(trial_data(1).(fn_meta{f}))
            avg_data(i).(fn_meta{f}) = unique([trial_data(cond_idx{i}).(fn_meta{f})]);
        else
            if size(trial_data(1).(fn_meta{f}),2) > 1
                u = unique(cat(1,trial_data(cond_idx{i}).(fn_meta{f})),'rows');
            else
                u = unique([trial_data(cond_idx{i}).(fn_meta{f})]);
            end
            avg_data(i).(fn_meta{f}) = u;
        end
    end
    
    % now loop along time signals to average
    for v = 1:length(time_vars)
        avg_data(i).(time_vars{v}) = mean(cat(3,trial_data(cond_idx{i}).(time_vars{v})),3);
        avg_data(i).([time_vars{v} '_std']) = std(cat(3,trial_data(cond_idx{i}).(time_vars{v})),[],3);
    end
    if avg_flag
        avg_data(i).is_average = true;
    end
    % add idx fields
    fn_idx = getTDfields(trial_data,'idx');
    
        for f = 1:length(fn_idx)
            if length(unique([trial_data(cond_idx{i}).(fn_idx{f})])) == 1
                avg_data(i).(fn_idx{f}) = trial_data(cond_idx{i}(1)).(fn_idx{f});
            end
        end

end

% restore logical order
avg_data = reorderTDfields(avg_data);


