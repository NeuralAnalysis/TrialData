%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getNorm(trial_data, params)
%
%   Takes the norm at each time point. Instead of params struct, can pass
% in signals.
%
% INPUTS:
%   trial_data : THE STRUCT
%   params     : parameter struct
%       .signals   : (string/cell) which signal(s) to use
%       .norm_name : (string/cell) field name to use for the stored norm.
%                      Default (empty) will just append '_norm'
%                      after the signal of interest. NOTE that if multiple
%                      signals are requested, this should be a cell array
%                      of names for each!
%                      
%   NOTE : can pass signals in as second input instead of params struct
%
% OUTPUTS:
%   trial_data : the struct (with added norm field)
%
% EXAMPLES:
%   To compute speed from velocity signals ('vel')
%       trial_data = getNorm(trial_data,struct('signals','vel','norm_name','speed');
%
% Written by Matt Perich. Updated October 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data  = getNorm(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals    =  ''; % which signals to use
norm_name  =  {''}; % the name of the normalized signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    error('Must define signals!');
else
    if ~isstruct(params)
        signals = params;
    else
        assignParams(who,params);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

signals = check_signals(trial_data,signals);

% must be a cell
if ~iscell(norm_name), norm_name = {norm_name}; end

% define all of the names to use
if isempty(norm_name) ||  ( length(norm_name) == 1 && isempty(norm_name{1}) )
    all_norm_names = signals;
    for iSig = 1:size(signals,1)
        all_norm_names{iSig} = [all_norm_names{iSig} '_norm'];
    end
end

%  check to make sure we have all the names
if length(norm_name) ~= size(signals,1)
    error('Not enough names for the fields! Need one per signal.');
end

% take norm for all trials and signals
for trial = 1:length(trial_data)
    for iSig = 1:size(signals,1)
        sig = trial_data(trial).(signals{iSig,1});
        n_time = size(sig,1);
        new_sig = zeros(n_time,1);
        for t = 1:n_time
            new_sig(t) = norm(sig(t,signals{iSig,2}));
        end
        
        trial_data(trial).(norm_name{iSig}) = new_sig;
    end
end