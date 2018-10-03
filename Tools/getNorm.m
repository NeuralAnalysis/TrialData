%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getNorm(trial_data,params)
%
%   Takes the norm at each time point. Instead of params struct, can pass
% in signals.
%
% Written by Matt Perich. Updated October 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data  = getNorm(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norm_name = 'norm'; % appends this after the signal name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(params)
    signals = params;
else
    assignParams(who,params);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

signals = check_signals(trial_data,signals);

for trial = 1:length(trial_data)
    for iSig = 1:size(signals,1)
        sig = trial_data(trial).(signals{iSig,1});
        n_time = size(sig,1);
        new_sig = zeros(n_time,1);
        for t = 1:n_time
            new_sig(t) = norm(sig(t,signals{iSig,2}));
        end
        
        trial_data(trial).([signals{iSig,1} '_' norm_name]) = new_sig;
    end
end