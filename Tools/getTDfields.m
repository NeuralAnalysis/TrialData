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
    case 'time'
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
    case 'cont' % same as 'time' but I exclude spikes
        cont_vars = fn(ismember(fn,cont_vars));
        t = size(trial_data(1).(cont_vars{1}),1);
        idx = false(1,length(fn));
        for ifn = 1:length(fn)
            idx(ifn) = size(trial_data(1).(fn{ifn}),1)==t;
        end
        % haven't thought of a cleverer way yet. I just hard coded a check
        % to exclude the types of neural things
        idx = idx & ...
            cellfun(@(x) isempty(x),strfind(fieldnames(trial_data),'_spikes'))' & ...
            cellfun(@(x) isempty(x),strfind(fieldnames(trial_data),'_pca'))' & ...
            cellfun(@(x) isempty(x),strfind(fieldnames(trial_data),'_gpfa'))' & ...
            cellfun(@(x) isempty(x),strfind(fieldnames(trial_data),'potent'))' & ...
            cellfun(@(x) isempty(x),strfind(fieldnames(trial_data),'null'))';
        fn = fn(idx);
    case 'spikes'
        fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'_spikes')));
    case 'arrays' % same as spikes but I only return the array name
        fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'_spikes')));
        fn = strrep(fn,'_spikes','')';
    case 'idx'
        fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'idx_')));
    otherwise
        warning('Category not recognized.')
        fn = {};
end