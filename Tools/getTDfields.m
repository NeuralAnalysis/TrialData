%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fn = getTDfields(trial_data,which_type)
%
%   Will return a cell array which lists all of the field names of a given
% type. Useful to list all spiking variables, kinematic variables, etc.
%
% INPUTS:
%   trial_data : the struct
%   which_type : (string) the type of field. Options:
%                   1) 'cont'        : continuous (pos, vel, force, etc)
%                   2) 'spikes'      : neural data fields
%                   3) 'arrays'      : similar to spikes but returns array names
%                   4) 'time'        : names of all time varying fields
%                   5) 'idx'         : name of all time index fields
%                   6) 'neural'      : any neural signals (e.g. M1_WHATEVER)
%                   7) 'unit_guides' : all unit_guides
%                   8) 'meta'        : all fields that are not time-varying or idx_
%
% OUTPUTS:
%   fn : the fieldnames of which_type
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn = getTDfields(trial_data,which_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont_vars = {'pos','vel','speed','acc','force','emg'}; % hard coded list of options
% these vars are common and known to be meta. Useful for edge case outlined below in time
meta_vars = {'trial_id','target_direction','target_center','bin_size','perturbation_info'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn = fieldnames(trial_data);
switch lower(which_type)
    case 'time' % all fields that are time-varying
        % find any signal that is KNOWN to be time-varying (defined above)
        % then find all signals that have the same number of rows
        % kinda hack-y but it works
        %   note: assumes rows are time and columns are variables
        cont_vars_here = fn(ismember(fn,cont_vars));
        % use the max over all trials so we have the lowest chance of
        % getting zero or one
        [t,trial_idx] = max(cellfun(@(x) size(x,1),{trial_data.(cont_vars_here{1})}));
        if t == 0
            error('Time variables have zero bins.');
        elseif t > 1
            idx = false(length(fn),1);
            for ifn = 1:length(fn)
                idx(ifn) = size(trial_data(trial_idx).(fn{ifn}),1)==t;
            end
            % in RW, the target center field has rows equal to the number
            % of targets, so hard code in an exclusion here
            idx = idx & ~ismember(fn,'target_center');
            fn = fn(idx);
        elseif t == 1
            % there is an edge case where if there is only one time bin it
            % returns everything. So, try and eliminate things through
            % process of elimination
            warning('Only one time bin for time-varying signals. getTDfields may be unreliable.');
            fn_idx = getTDfields(trial_data,'idx');
            fn_ug  = getTDfields(trial_data,'unit_guides');
            fn_spikes = getTDfields(trial_data,'spikes');
            bad_idx = zeros(length(fn),1);
            for i = 1:length(fn)
                if ischar(trial_data(1).(fn{i}))
                    % character fields are not time-varying
                    bad_idx(i) = 1;
                elseif size(trial_data(1).(fn{i}),2) == 1
                    % probably not time varying if it's only one point
                    bad_idx(i) = 1;
                elseif size(trial_data(1).(fn{i}),1) ~= t
                    bad_idx(i) = 1;
                end
            end
            bad_idx = bad_idx | ismember(fn,fn_idx) | ismember(fn,fn_ug) | ismember(fn,meta_vars);
            good_idx = zeros(length(fn),1);
            for i = 1:length(cont_vars)
                good_idx = good_idx | cellfun(@(x) ~isempty(x),strfind(fn,cont_vars{i}));
            end
            good_idx = good_idx | ismember(fn,fn_spikes);
            fn = fn(good_idx | ~bad_idx);
        end
        
    case 'meta'% all fields that are NOT time-varying or idx_
        % find any signal that is KNOWN to be time-varying (defined above)
        % then find all signals that have the same number of rows
        % kinda hack-y but it works
        %   note: assumes rows are time and columns are variables
        fn_time = getTDfields(trial_data,'time');
        fn_idx  = getTDfields(trial_data,'idx');
        fn_ug   = getTDfields(trial_data,'unit_guides');
        idx = ismember(fn,fn_time) | ismember(fn,fn_idx) | ismember(fn,fn_ug);
        fn = fn(~idx);
    case 'cont' % same as 'time' but I exclude neural
        fn_time = getTDfields(trial_data,'time');
        fn_neural = getTDfields(trial_data,'neural');
        fn_unit_guides = getTDfields(trial_data,'unit_guides');
        fn = fn(ismember(fn,fn_time) & ~ismember(fn,fn_neural) & ~ismember(fn,fn_unit_guides));
    case 'spikes' % just the _spikes fields
        fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'_spikes')));
    case 'unit_guides'
        fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'_unit_guide')));
    case 'arrays' % same as spikes but I only return the array name
        fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'_spikes')));
        fn = strrep(fn,'_spikes','')';
    case 'neural' % anything that is neural derived (e.g. M1_spikes and M1_pca)
        arrays = getTDfields(trial_data,'arrays');
        fn = fieldnames(trial_data);
        neural_idx = zeros(length(fn),1);
        for array = 1:length(arrays)
            neural_idx = neural_idx | ...
                (cellfun(@(x) ~isempty(regexp(x,[arrays{array} '_'],'ONCE')),fn) & ...
                ~cellfun(@(x) ~isempty(regexp(x,[arrays{array} '_unit_guide'],'ONCE')),fn));
        end
        fn = fn(neural_idx);
    case 'idx' % any idx_ field
        fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'idx_')));
    otherwise
        disp('Category not recognized.')
        fn = {};
end