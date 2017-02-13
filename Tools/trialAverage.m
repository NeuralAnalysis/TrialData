%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [avg_data, cond_idx] = trialAverage(trial_data, conditions, params)
%
%   Returns struct where each entry is an average for all time-varying
% signals across trials for some condition. This function can interpolate
% to stretch/shrink all trials to the same number of points, but it does
% not trim data. This can be easily done with truncateAndBin.
%
% INPUTS:
%   trial_data : the struct
%   conditions : (string or cell array) the field name(s) in trial_data
%                   avg_data will have an entry for each unique combo
%   params     : struct with parameters
%     .do_stretch : (bool) whether to stretch/shrink trials to uniform length
%                       if false, all trials must have same number of points
%     .num_samp   : (int) how many time points to use for interpolation
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
if nargin > 2
    eval(structvars(length(fieldnames(params)),params)); %overwrite parameters
end
if ~iscell(conditions), conditions = {conditions}; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get list of time-varying signals that we will average over
time_vars = getTDfields(trial_data,'time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time warp each trial to the same number of points, if desired
if do_stretch
    for trial = 1:length(trial_data)
        for iVar = 1:length(time_vars)
            temp = trial_data(trial).(time_vars{iVar});
            trial_data(trial).(time_vars{iVar}) = interp1(1:size(temp,1),temp,linspace(1,size(temp,1),num_samp));
        end
    end
end

if length(unique(cellfun(@(x) size(x,1),{trial_data.pos}))) ~= 1
    error('Trials are not uniform length. Do this with time stretching option or truncateAndBin.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop along conditions and get unique values for each
cond_vals = cell(1,length(conditions));
for iCond = 1:length(conditions)
    % unique doesn't work on numeric cell arrays for some reason
    if ischar(trial_data(1).(conditions{iCond}))
        uc = unique({trial_data.(conditions{iCond})});
    elseif isnumeric(trial_data(1).(conditions{iCond}))
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
cond_idx = cell(1,num_conds);
avg_data = repmat(struct(),1,num_conds);
for i = 1:num_conds
    avg_data(i).monkey = trial_data(1).monkey;
    avg_data(i).date = trial_data(1).date;
    avg_data(i).task = trial_data(1).task;
    avg_data(i).bin_size = trial_data(1).bin_size;
    
    func_in = cell(1,2*size(all_conds,2));
    for iCond = 1:size(all_conds,2)
        func_in{2*(iCond-1)+1}   = conditions{iCond};
        func_in{2*(iCond-1)+2} = cond_vals{iCond}{all_conds(i,iCond)};
        avg_data(i).(conditions{iCond}) = cond_vals{iCond}{all_conds(i,iCond)};
    end
    cond_idx{i}=getTDidx(trial_data,func_in);
    
    % now loop along time signals to average
    for v = 1:length(time_vars)
        avg_data(i).(time_vars{v}) = mean(cat(3,trial_data(cond_idx{i}).(time_vars{v})),3);
    end
    
end


