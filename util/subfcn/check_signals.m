%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function signals = check_signals(trial_data,signals)
%
%   Basic utility function to parse out the "signals" parameters typically
% used for getPCA, getModel, etc, and return something that the function is
% sure to accept.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signals = check_signals(trial_data,signals)
% make sure it's a cell (you can pass in a string)
if ~iscell(signals), signals = {signals}; end
% see if it's 1D or 2D
if any(size(signals)==1) % it may be 1-D name list
    if numel(signals) == 1 % it's only a single name, so use all columns
        signal_idx = {1:size(trial_data(1).(signals{1}),2)};
    elseif numel(signals) == 2 && (~ischar(signals{2}) || strcmpi(signals{2},'all')) % second entry is an idx
        signal_idx = signals(2);
        signals = signals(1);
    elseif numel(signals) > 2 % it's a list of names
        signal_idx = repmat({'all'},length(signals),1);
        if size(signals,2) == 1, signals = signals'; end
    end
    % return in the expected format
    signals = [signals, signal_idx];
elseif size(signals,2) > 2
    error('Signals must be cell vector of names or Nx2 cell array of name/idx pairs.');
end %otherwise, things look good!

% Check for 'all'
% figure out how many signals there will be
for i = 1:size(signals,1)
    % if second entry is 'all', use all
    if ischar(signals{i,2})
        if strcmpi(signals{i,2},'all') % replace with all indices
            signals{i,2} = 1:size(trial_data.(signals{i,1}),2);
        else, error('Not sure what this signal idx is...');
        end
    end
end