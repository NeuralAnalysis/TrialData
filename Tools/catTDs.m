%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = catTDs(varargin)
%
%   Concatenates the TD structs, even if they don't have the same fields
%
% INPUTS:
%   varargin : any number of trial_data structs to be concatenated
%
% EXAMPLE:
%   all_the_tds = catTDs(td1,td2,td3);
%
% Written by Matt Perich. Updated August 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function master_td = catTDs(varargin)
% do some checks of varargin (e.g. are they structs)
for iTD = 1:length(varargin)
    if ~isstruct(varargin{iTD})
        error('Must provide only trial_data structs as input');
    end
end

% start the concatenation
if length(varargin) > 1
    % initialize the master
    master_td = [];
    for iTD = 1:length(varargin)
        
        if isempty(master_td)
            master_td = varargin{iTD};
        end
        
        if ~isempty(master_td)
            td = varargin{iTD};
            
            fn_master = fieldnames(master_td);
            fn_new = fieldnames(td);
            
            % see what is missing from the master, and populate it with NaN
            fn_miss = setdiff(fn_new,fn_master);
            for i = 1:length(fn_miss)
                for trial = 1:length(master_td)
                    master_td(trial).(fn_miss{i}) = NaN;
                end
            end
            
            % see what is missing from the new, and populate it with NaN
            fn_miss = setdiff(fn_master,fn_new);
            for i = 1:length(fn_miss)
                for trial = 1:length(td)
                    td(trial).(fn_miss{i}) = NaN;
                end
            end
            
            master_td = [master_td,td];
        end
    end
else
    disp('Only one trial_data struct. Nothing to concatenate...');
    master_td = varargin{1};
end
