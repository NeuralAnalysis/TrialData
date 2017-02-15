%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fn = getTDfields(trial_data,which_type)
%
%   Will return a cell array which lists all of the field names of a given
% type. Useful to list all spiking variables, kinematic variables, etc.
%
% INPUTS:
%   trial_data : the struct
%   which_type : (string) the type of field. Options:
%                   1) 'cont'   : continuous (pos, vel, force, etc)
%                   2) 'spikes' : neural data fields
%                   3) 'arrays' : similar to spikes but returns array names
%                   4) 'time'   : names of all time varying fields
%                   5) 'idx'    : name of all time index fields
%                   6) 'neural'  : any neural signals (e.g. M1_WHATEVER)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn = fieldnames(trial_data);
switch lower(which_type)
    case 'time' % all fields that are time-varying
        % find any signal that is KNOWN to be time-varying (defined above)
        % then find all signals that have the same number of rows
        % kinda hack-y but it works
        %   note: assumes rows are time and columns are variables
        cont_vars = fn(ismember(fn,cont_vars));
        t = size(trial_data(1).(cont_vars{1}),1);
        idx = false(1,length(fn));
        for ifn = 1:length(fn)
            idx(ifn) = size(trial_data(1).(fn{ifn}),1)==t;
        end
        fn = fn(idx);
    case 'cont' % same as 'time' but I exclude neural
        cont_vars = fn(ismember(fn,cont_vars));
        t = size(trial_data(1).(cont_vars{1}),1);
        idx = false(length(fn),1);
        for ifn = 1:length(fn)
            idx(ifn) = size(trial_data(1).(fn{ifn}),1)==t;
        end
        fn_neural = getTDfields(trial_data,'neural');
        fn = fn(idx & ~ismember(fn,fn_neural));
    case 'spikes' % just the _spikes fields
        fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'_spikes')));
    case 'arrays' % same as spikes but I only return the array name
        fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'_spikes')));
        fn = strrep(fn,'_spikes','')';
    case 'neural' % anything that is neural derived (e.g. M1_spikes and M1_pca)
        arrays = getTDfields(trial_data,'arrays');
        fn = fieldnames(trial_data);
        neural_idx = zeros(length(fn),1);
        for array = 1:length(arrays)
            neural_idx = neural_idx | cellfun(@(x) ~isempty(regexp(x,'_spikes','ONCE')),fn);
        end
        fn = fn(neural_idx);
    case 'idx' % any idx_ field
        fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'idx_')));
    otherwise
        disp('Category not recognized.')
        fn = {};
end