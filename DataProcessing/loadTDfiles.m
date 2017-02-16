%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function master_td = loadTDfiles(varargin)
%
%   Loads an arbitrary number of trial_data files from disk and appends
% them into one structure. If any fields don't match, will fill in the
% missing bits with NaN entries.
%
% INPUTS:
%   filenames : a cell array of filenames
%   varargin : any number of function call inputs
%               Format 1: {'functionName',params}
%                       Note: the params struct should be made as you
%                             would for the real function call.
%               Format 2: 'functionName', will assume no params
%               Format 3: {'functionName',var1,var2,etc} for arbitrary vars
%
% OUTPUTS:
%   master_td : the master trial_data struct with all files appended and
%               all functions in varargin executed
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function master_td = loadTDfiles(filenames,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make sure the inputs make sense
if ~iscell(filenames), filenames = {filenames}; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(filenames) > 1
    disp(['Loading file ' num2str(1) ' of ' num2str(length(filenames)) '.']);
    if exist(filenames{1},'file')
        load(filenames{1});
    else
        error([filenames{1} ' not found.']);
    end
    
    if nargin > 1
        for iFun = 1:length(varargin)
            trial_data = run_func(trial_data,varargin{iFun});
        end
    end
    
    master_td = trial_data;
    
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
                trial_data = run_func(trial_data,varargin{iFun});
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
    load(filenames{1});
    master_td = trial_data;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = run_func(trial_data,funcCall)
if ~iscell(funcCall), funcCall = {funcCall}; end
if exist(funcCall{1},'file') || exist(funcCall{1},'builtin')
    % build function call
    f = eval(['@' funcCall{1}]);
    % call it
    if length(funcCall) > 1
        if isstruct(funcCall{2})
            trial_data = f(trial_data,funcCall{2});
        else
            input_str = [];
            for i = 2:length(funcCall)
                if isnumeric(funcCall{i})
                    in_var = num2str(funcCall{i});
                else
                    in_var = funcCall{i};
                end
                input_str = [input_str ',' in_var];
            end
            eval(['trial_data = f(trial_data' input_str ');']);
        end
    else
        trial_data = f(trial_data);
    end
else
    warning([funcCall ' function not found. Skipping...']);
end

% restore logical order
trial_data = reorderTDfields(trial_data);
end