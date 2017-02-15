%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho,sort_idx] = pairwiseCorr(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [rho,sort_idx] = pairwiseCorr(trial_data,params)
%
%   Compute pairwise correlations between neurons in _spikes.
%   Can order them to show structure, a la Cunningham.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%     .arrays         : which units (can be cell array with multiple)
%     .trial_idx      : which trials to use (default: all)
%     .neurons        : which neurons to use (default: all) Note: for multiple arrays,
%                       neurons should be cell array with indices for each array
%     .cluster_order  : flag to reorder cells to show structure (default: false)
%     .array_cluster  : extra flag to keep arrays groups separately when reordering
%     .do_norm        : normalize each row/col to go from -1 to 1
% 
% OUTPUTS:
%   rho        : pairwise correlation matrix for _spikes
%   sort_idx   : indices of original _spikes inputs for sort
%
% Written by Matt Perich. Updated Feb 2017.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS
arrays          =  [];
trial_idx       =  1:length(trial_data);
neurons         =  [];
cluster_order   =  false;
array_cluster   =  false;
do_norm         =  false;
if nargin > 1, assignParams(who,params); end % overwrite defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if ~iscell(arrays), arrays = {arrays}; end
if isempty(neurons)
    neurons = cell(1,length(arrays));
    for i = 1:length(arrays), neurons{i} = 1:size(trial_data(1).([arrays{i} '_spikes']),2); end
end
if ~iscell(neurons), neurons = {neurons}; end
if isempty(arrays), error('Must give array'); end
trial_data = trial_data(trial_idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build master matrix of spiking
fr = [];
for array = 1:length(arrays)
    temp = cat(1,trial_data.([arrays{array} '_spikes']));
    fr = [fr, temp(:,neurons{array})];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get pairwise correlations, and replace diagonal with zeros
rho = corr(fr).*(-1*eye(size(fr,2))+ones(size(fr,2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize each row/column
if do_norm
    [temp1,temp2] = deal(zeros(size(rho)));
    temp1(rho > 0) = rho(rho > 0);
    temp1 = temp1./repmat(max(rho),size(rho,1),1);
    
    temp2(rho < 0) = rho(rho < 0);
    temp2 = temp2./repmat(min(rho),size(rho,1),1);
    rho = temp1 + temp2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reorder to highlight structure
if cluster_order
    if array_cluster
        temp = cell(1,length(arrays));
        for array = 1:length(arrays)
            temp{array} = cluster_rho(rho,1:neurons{array}(end));
        end
        
    else
        sort_idx = cluster_rho(rho);
    end
    rho = rho(sort_idx,sort_idx);
else
    sort_idx = 1:size(rho,2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sort_idx = cluster_rho(rho)

% z = squareform(pdist(rho));
c=cluster(linkage(rho),'cutoff',1,'depth',3);
u = unique(c);

sort_idx = [];
for i = 1:length(u)
    idx = find(c==u(i));
    sort_idx = [sort_idx; idx];
end

end
