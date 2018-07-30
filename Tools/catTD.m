function master_td = catTD(varargin)
% concatenates the TD structs, even if they don't have the same fields

% do some checks of varargin (e.g. are they structs)

% initialize the master
master_td = [];

for iS = 1:length(varargin)
    
    if isempty(master_td)
        master_td = varargin{iS};
    end
    
    if ~isempty(master_td)
        td = varargin{iS};
        
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
