%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = stretchSignals(trial_data,params)
%
%   Time warps signals to be the same length. Useful for averaging, etc.
%   Note: for idx_ fields, adjusts to new sampling length and uses ceil to
%   make it an integer. So the precision of these idx_ fields now depends
%   on the num_samp value chosen.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = stretchSignals(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
num_samp    =  100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end
% get list of time-varying signals that we will average over
time_vars = getTDfields(trial_data,'time');
idx_vars = getTDfields(trial_data,'idx');

for trial = 1:length(trial_data)
        
    for iVar = 1:length(idx_vars)
        og_size = size(trial_data(trial).(time_vars{1}),1);
        trial_data(trial).(idx_vars{iVar}) = ceil(num_samp*trial_data(trial).(idx_vars{iVar})/og_size);
    end
    
    for iVar = 1:length(time_vars)
        temp = trial_data(trial).(time_vars{iVar});
        trial_data(trial).(time_vars{iVar}) = interp1(1:size(temp,1),temp,linspace(1,size(temp,1),num_samp));
    end
end
