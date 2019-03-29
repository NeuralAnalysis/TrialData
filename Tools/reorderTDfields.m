%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = reorderTDfields(trial_data)
%
% This is my total OCD-appeasing function. It will:
%   1) Put metadata at the top
%   2) Put idx_ fields second, in order of value for each trial
%   3) put continuous/time signals third (EMG last)
%   4) put spikes and neural fields last (and group by array name)
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = reorderTDfields(trial_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
% Get names for the different types
fn        = fieldnames(trial_data);
fn_idx    = getTDfields(trial_data,'idx');
fn_neural = [getTDfields(trial_data,'neural'); getTDfields(trial_data,'unit_guides')];
fn_cont   = [getTDfields(trial_data,'cont'); getTDfields(trial_data,'labels')];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do some checks for bad field names
if ismember('all',fn)
    warning('TD struct has a field named ''all''... this could break some things, and is not a very informative field anyway. Consider changing it.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get indices for everything

% the easy ones
meta_idx = find(~(ismember(fn,fn_idx) | ismember(fn,fn_neural) | ismember(fn,fn_cont)));
neural_idx = find(ismember(fn,fn_neural));

% sort idx indices by increasing bin count
idx_idx = find(ismember(fn,fn_idx));
% DO SORT
v = zeros(1,length(idx_idx));
for i = 1:length(idx_idx)
    temp = trial_data(1).(fn{idx_idx(i)});
    if isempty(temp) || any(isnan(temp))
        v(i) = Inf;
    else
        v(i) = temp(1); % in case there are multiple
    end
end
[~,I] = sort(v);
idx_idx = idx_idx(I);


%%%% THIS CODE IS VERY VERY INEFFICIENT BUT IT KINDA WORKS. I JUST WROTE IT
%%%% AS QUICKLY AS POSSIBLE TO MAKE IT WORK. SOMEDAY, SOMEONE SHOULD MAKE
%%%% IT MORE CLEVER AND CLEAN
% get continuous
cont_idx = find(ismember(fn,fn_cont));
% we want to reorder such that _shift fields are near their originals
shift_idx_all = fn_cont(~cellfun(@isempty,regexp(fn_cont,'_shift','ONCE')));
if ~isempty(shift_idx_all)
    fn_cont_new = fn_cont;
    for i = 1:length(shift_idx_all)
        reorder_idx = (1:length(fn_cont_new))';
        
        % get the current index of this _shift field
        shift_idx = find(~cellfun(@isempty,regexp(fn_cont_new,shift_idx_all{i},'ONCE')));
        
        the_val = shift_idx;
        og_name = fn_cont_new{the_val};
        new_name = strsplit(og_name,'_shift'); % effectively removes the _shift
        new_name = cat(2,new_name{~cellfun(@isempty,new_name)});
        
        % find what the original name was and get all possible ones
        og_idx = setdiff(find(~cellfun(@isempty,regexp(fn_cont_new,new_name,'ONCE'))),shift_idx);
        if isempty(og_idx)
            error('Difficulty sorting out the _shift fields. No idea why!');
        end
        og_idx = og_idx(end); % get the last one
        
        reorder_idx(shift_idx) = NaN;
        reorder_idx = [reorder_idx(1:og_idx); the_val; reorder_idx(og_idx+1:end)];
        reorder_idx(isnan(reorder_idx)) = [];
        
        fn_cont_new = fn_cont_new(reorder_idx);
    end
    % now map this new order onto the master fieldname list
    cont_idx_new = zeros(size(cont_idx));
    for i = 1:length(fn_cont)
        temp_idx = strcmpi(fn_cont,fn_cont_new{i});
        cont_idx_new(i) = cont_idx(temp_idx);
    end
    cont_idx = cont_idx_new;
end

% put EMG at end of continuous
%   Because, aesthetically, I like having emg_names after the rest
if isfield(trial_data,'emg')
    emg_idx = cellfun(@(x) ~isempty(x),strfind(fn_cont,'emg'));
    cont_idx = [cont_idx; cont_idx(emg_idx)];
    cont_idx(emg_idx) = [];
end

% do reordering
master_idx = [meta_idx; idx_idx; cont_idx; neural_idx];

% order them
trial_data = orderfields(trial_data,master_idx);
