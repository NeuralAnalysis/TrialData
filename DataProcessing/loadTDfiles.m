%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [master_td, params] = loadTDfiles(varargin)
%
%   Loads an arbitrary number of trial_data files from disk and appends
% them into one structure. If any fields don't match, will fill in the
% missing bits with NaN entries.
%
% INPUTS:
%   filenames : a cell array of filenames
%   varargin : any number of function call inputs
%               Format 1: {@functionName,params}
%                       Note: the params struct should be made as you
%                             would for the real function call.
%               Format 2: @functionName, will assume no params
%               Format 3: {@functionName,var1,var2,etc} for arbitrary vars
%
% OUTPUTS:
%   master_td : the master trial_data struct with all files appended and
%               all functions in varargin executed
%   params    : struct of info
%       .func_calls : the function calls (varargin)
%       .git_info   : the git details for the current TrialData repo
%       .extra_outs : extra outputs from functions that have them
%                       Note: this field will be missing if no functions
%                             return outputs, for the sake of cleanliness
%
% EXAMPLES:
%   e.g. to bin data and trim it
%   [td,params] = loadTDfiles(filename, {@binTD,5}, {@trimTD,{'idx_go_cue',0},{'idx_go_cue',30}});
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [master_td, params] = loadTDfiles(filenames,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make sure the inputs make sense
if ~iscell(filenames), filenames = {filenames}; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

extra_outs = cell(length(filenames),length(varargin));

% Load the first file to get it started. Will then loop along any extra
% files if desired
disp(['Loading file ' num2str(1) ' of ' num2str(length(filenames)) '.']);
[master_td,extra_outs(1,:)] = do_the_loading(filenames{1},varargin);

if length(filenames) > 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop along filenames
    for file = 2:length(filenames)
        disp(['Loading file ' num2str(file) ' of ' num2str(length(filenames)) '.']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [td,extra_outs(file,:)] = do_the_loading(filenames{file},varargin);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check fieldnames of new file against others
        master_fn = fieldnames(master_td);
        fn = fieldnames(td);
        [~,Ia,Ib] = setxor(master_fn,fn);
        % any fields not in common get filled in with NaNs
        if ~isempty(Ia)
            for i = 1:length(Ia)
                disp(['Field ''' master_fn{Ia(i)} ''' is missing from File ' num2str(file) '.']);
                [td.(master_fn{Ia(i)})] = deal(NaN);
            end
        end
        if ~isempty(Ib)
            for i = 1:length(Ib)
                disp(['File ' num2str(file) ' has new field ''' fn{Ib(i)} '''.']);
                [master_td.(fn{Ib(i)})] = deal(NaN);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % concatenate and move on
        master_td = [master_td, td];
    end
else
    disp('Only one filename provided. Returning one file.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restore logical order
master_td = reorderTDfields(master_td);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% package up params if desired
if nargout > 1
    params.filenames      = filenames;
    params.func_calls     = varargin;
    params.git_info       = getGitInfo();
    if ~all(all(cellfun(@isempty,extra_outs)))
        params.extra_outs = extra_outs;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function loads stuff
function [trial_data,extra_outs] = do_the_loading(filename,func_calls)
if ~ischar(filename), error('filename must be a string with a path to the trial_data file.'); end
if exist(filename,'file')
    load(filename);
else
    error([filename ' not found.']);
end
% run processing functions if provided
if nargin > 1
    extra_outs = cell(1,length(func_calls));
    for iFun = 1:length(func_calls)
        [trial_data,extra_outs{1,iFun}] = run_func(trial_data,func_calls{iFun});
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function runs each of the input handles and parameters
function [trial_data,out_params] = run_func(trial_data,funcCall)
if ~iscell(funcCall), funcCall = {funcCall}; end
out_params = [];

% first entry is function call!
f = funcCall{1};
fh = functions(f);
% make sure it exists
if ~isempty(fh.file)
    % check if there are parameters too
    if length(funcCall) > 1
        % there is a special case for getTDidx
        if strcmpi(fh.function,'getTDidx')
            [~,trial_data] = f(trial_data,funcCall{2:end});
        else
            if nargout(f) > 1 % second output is stored in params
                [trial_data,out_params] = f(trial_data,funcCall{2:end});
            else
                trial_data = f(trial_data,funcCall{2:end});
            end
        end
    else % no parameters
        trial_data = f(trial_data);
    end
else
    error([funcCall ' function not found...']);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
