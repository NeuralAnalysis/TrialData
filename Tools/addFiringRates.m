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
array = '';
method = 'averageInBin';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here are some parameters that you can overwrite that aren't documented
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end

switch lower(method)
    case 'averageinbin'
        for trial_idx = 1:length(trial_data)
            spikes = trial_data(trial_idx).([array '_spikes']);
        
            bin_size = trial_data(trial_idx).bin_size;
            FR = spikes/bin_size;
            trial_data(trial_idx).([array '_FR']) = FR;
        end
    otherwise
        error('That FR calculating method is not implemented')
end

td = trial_data;
