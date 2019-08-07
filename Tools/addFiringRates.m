%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function td = addFiringRates(trial_data,params)
% 
% Adds computes firing rates in Hz on given array, using given method
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .array          : which array to calculate FRs for
%       .method         : which method
%                           default: 'averageInBin'
%
% OUTPUTS:
%   td  : new TrialData struct with FRs added
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function td = addFiringRates(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
array         = '';
method        = 'averageInBin';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here are some parameters that you can overwrite that aren't documented
field_extra   =  {'_FR'};
verbose       = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if isstruct(params) % overwrite parameters
        assignParams(who,params);
    else % assumes you  passed in array
        array = params;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data  =  check_td_quality(trial_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check output field addition
field_extra  = check_field_extra(field_extra,array);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
switch lower(method)
    case 'averageinbin'
        for trial = 1:length(trial_data)
            spikes = trial_data(trial).([array '_spikes']);
        
            bin_size = trial_data(trial).bin_size;
            FR = spikes/bin_size;
            trial_data(trial).([array field_extra{1}]) = FR;
        end
    otherwise
        error('That FR calculating method is not implemented')
end

td = trial_data;
