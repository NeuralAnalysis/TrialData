%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data, signal_means] = centerTD(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This will center (ie mean-subtract) any (or all) time varying signals
% for the provided TD. Nothing fancy.
%
%   Can provide signals as second input instead of params struct.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .signals : which signal(s) to center
%
% OUTPUTS:
%   trial_data  : the struct with centered signals
%   signal_means : cell array of means vectors for all signals
%
% Written by Matt Perich. Updated March 2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals = '';
if nargin > 1
    if isstruct(params)
        assignParams(who,params);
    else
        signals = params;
    end
else
    signals = getTDfields(trial_data,'time');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals = check_signals(trial_data,signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

signal_means = cell(1,size(signals,1));
for iSig = 1:size(signals,1)
    m = mean(getSig(trial_data,signals(iSig,:)),1);
    for trial = 1:length(trial_data)
        temp = getSig(trial_data(trial),signals(iSig,:));
        trial_data(trial).(signals{iSig,1}) = temp - repmat(m,size(temp,1),1);
    end
    signal_means{iSig} = m;
end
