%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x = get_vars(trial_data, signals)
%
%   Boring utility subfunction used by things like getModel to turn a
% signals input into a matrix of data. The signals parameter is the
% standard {'NAME',idx; etc} format used by getPCA and getModel, so this
% could be a generally useful function to you...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = get_vars(trial_data, signals)
% figure out how many signals there will be
idx = cell(1,size(signals,1));
for i = 1:size(signals,1)
    idx{i} = signals{i,2};
end

% preallocate, because with massively large datasets its worth it
fn = getTDfields(trial_data,'time');
x = zeros(size(cat(1,trial_data.(fn{1})),1),sum(cellfun(@(x) length(x),idx)));
count = 0;
% piece everything together
for i = 1:size(signals,1)
    temp = cat(1,trial_data.(signals{i,1}));
    x(:,count+(1:length(idx{i}))) = temp(:,idx{i});
    count = count + length(idx{i});
end