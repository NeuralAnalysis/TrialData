%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = stretchSignals(trial_data,params)
%
%   Time warps signals to be the same length. Useful for averaging, etc.
%   Note: for idx_ fields, adjusts to new sampling length and uses ceil to
%   make it an integer. So the precision of these idx_ fields now depends
%   on the 'samples' value chosen.
%
% INPUTS:
%   trial_data : the struct
%   params     : params struct
%       .samples  : how many samples to use in interpolation (default: 100)
% OUTPUTS:
%   trial_data : the struct with signals that are time-stretched
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = stretchSignals(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
samples    =  100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented extra parameters
verbose      =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if ~isstruct(params), error('Second input must be params struct'); end
    if isfield(params,'num_samp')
        warning('stretchSignals: ''num_samp'' input has been changed to ''samples'' for ease of use. Please change your code. For now the parameter is  overwritten but this warning and the overwriting will be removed in a future release');
        samples = params.num_samp;
    end
    assignParams(who,params);
end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
% get list of time-varying signals that we will average over
time_vars = getTDfields(trial_data,'time');
idx_vars = getTDfields(trial_data,'idx');

for trial = 1:length(trial_data)
        
    for iVar = 1:length(idx_vars)
        og_size = size(trial_data(trial).(time_vars{1}),1);
        trial_data(trial).(idx_vars{iVar}) = ceil(samples*trial_data(trial).(idx_vars{iVar})/og_size);
    end
    
    for iVar = 1:length(time_vars)
        temp = trial_data(trial).(time_vars{iVar});
        temp = interp1(1:size(temp,1),temp,linspace(1,size(temp,1),samples));
        if size(temp,1) == 1 && size(temp,2) ~= 1
            temp = temp';
        end
        trial_data(trial).(time_vars{iVar}) = temp;
    end
end
