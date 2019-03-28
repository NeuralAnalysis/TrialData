%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x = get_vars(trial_data, signals)
%
%   Boring utility subfunction used by things like getModel to turn a
% signals input into a matrix of data. The signals parameter is the
% standard {'NAME',idx; etc} format used by e.g. dimReduce and getModel, so
% this could be a generally useful function to you...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = get_vars(trial_data, signals)
% figure out how many signals there will be
idx = cell(1,size(signals,1));
for iSig = 1:size(signals,1)
    idx{iSig} = signals{iSig,2};
end

% preallocate, because with massively large datasets its worth it
x = zeros(size(cat(1,trial_data.(signals{1,1})),1),sum(cellfun(@(x) length(x),idx)));
count = 0;
% piece everything together
for iSig = 1:size(signals,1)
    temp = cat(1,trial_data.(signals{iSig,1}));
    x(:,count+(1:length(idx{iSig}))) = temp(:,idx{iSig});
    count = count + length(idx{iSig});
end