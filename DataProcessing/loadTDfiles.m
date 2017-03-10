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
%       .git_hash   : the git hash for the current TrialData repo
%       .extra_outs : extra outputs from functions that have them
%                       Note: this field will be missing if no functions
%                       return outputs, for the sake of cleanliness
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


disp(['Loading file ' num2str(1) ' of ' num2str(length(filenames)) '.']);
if exist(filenames{1},'file')
    load(filenames{1});
else
    error([filenames{1} ' not found.']);
end
extra_outs = cell(length(filenames),length(varargin));

if nargin > 1
    for iFun = 1:length(varargin)
        [trial_data,extra_outs{1,iFun}] = run_func(trial_data,varargin{iFun});
    end
end

master_td = trial_data;

if length(filenames) > 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop along filenames
    for file = 2:length(filenames)
        disp(['Loading file ' num2str(file) ' of ' num2str(length(filenames)) '.']);
        
        if exist(filenames{file},'file')
            load(filenames{file});
        else
            error([filenames{file} ' not found.']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Now run the arbitrary functions
        if nargin > 1
            for iFun = 1:length(varargin)
                [trial_data,extra_outs{file,iFun}] = run_func(trial_data,varargin{iFun});
            end
        end
        
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % concatenate and move on
        master_td = [master_td, trial_data];
    end
else
    disp('Only one filename provided. Returning one file.');
end

% restore logical order
master_td = reorderTDfields(master_td);

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
        if isstruct(funcCall{2}) % it's a params struct
            if nargout(f) > 1
                [trial_data,out_params] = f(trial_data,funcCall{2});
            else
                trial_data = f(trial_data,funcCall{2});
            end
        else % it's a series of inputs
            input_str = [];
            % loop along all of the inputs and build the function call
            for i = 2:length(funcCall)
                if isnumeric(funcCall{i}) % is it a number?
                    in_var = num2str(funcCall{i});
                else %no...
                    in_var = funcCall{i};
                    if ischar(in_var) % is it a string?
                        in_var = ['''' in_var ''''];
                    elseif iscell(in_var) % is it a cell? Must break it out
                        temp = '{';
                        for j = 1:length(in_var)
                            if ischar(in_var{j})
                                temp = [temp '''' in_var{j} ''''];
                            elseif isnumeric(in_var{j})
                                temp = [temp num2str(in_var{j})];
                            end
                            if j < length(in_var), temp = [temp ',']; end
                        end
                        in_var = [temp '}'];
                    end
                end
                input_str = [input_str ',' in_var];
            end
            % there is a special case for getTDidx
            if strcmpi(fh.function,'getTDidx')
                eval(['[~,trial_data] = f(trial_data' input_str ');']);
            elseif nargout(f) > 1
                eval(['[trial_data,out_params] = f(trial_data' input_str ');']);
            else
                eval(['trial_data = f(trial_data' input_str ');']);
            end
        end
    else % easy if there are no parameters
        trial_data = f(trial_data);
    end
else
    error([funcCall ' function not found...']);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
