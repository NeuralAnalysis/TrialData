%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function master_td = loadTDfiles(varargin)
%
%   Loads an arbitrary number of trial_data files from disk and appends
% them into one structure. If any fields don't match, will fill in the
% missing bits with NaN entries.
%
% INPUTS:
%   varargin : either a cell array of filenames or any number of strings
%
% OUTPUTS:
%   master_td : the master trial_data struct with all files appended
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function master_td = loadTDfiles(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make sure the inputs make sense
if length(varargin) == 1 && iscell(varargin{1}), varargin = varargin{1}; end
if ~iscell(varargin), varargin = {varargin}; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(varargin) > 1
    disp(['Loading file ' num2str(1) ' of ' num2str(length(varargin)) '.']);
    load(varargin{1});
    master_td = trial_data;
    for file = 2:length(varargin)
        disp(['Loading file ' num2str(file) ' of ' num2str(length(varargin)) '.']);
        load(varargin{file})

        % check fieldnames of new file against others
        master_fn = fieldnames(master_td);
        fn = fieldnames(trial_data);
        [~,Ia,Ib] = setxor(master_fn,fn);
        % any fields not in common get filled in with NaNs
        if ~isempty(Ia)
            for i = 1:length(Ia)
                disp(['Field ''' master_fn{Ia(i)} ''' is missing from File ' num2str(file) '.']);
                [trial_data.(master_fn{Ia(i)})] = deal(NaN);
            end
        end
        if ~isempty(Ib)
            for i = 1:length(Ib)
                disp(['File ' num2str(file) ' has new field ''' fn{Ib(i)} '''.']);
                [master_td.(fn{Ib(i)})] = deal(NaN);
            end
        end
        
        % concatenate and move on
        master_td = [master_td, trial_data];
    end
else
    disp('Only one filename provided. Returning one file.');
    load(varargin{1});
    master_td = trial_data;
end