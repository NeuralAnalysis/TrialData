%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [rho,sort_idx] = pairwiseCorr(trial_data,params)
%
%   Compute pairwise correlations between neurons in _spikes.
%   Can order them to show structure, a la Cunningham.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%     .signals        : which signals (format as in getPCA or getModel)
%     .method         : which method (default @corr, could do @cov for ex)
%     .trial_idx      : which trials to use (default: all)
%     .cluster_order  : flag to reorder cells to show structure (default: false)
%     .cluster_arrays : flag treat arrays differently when clustering
%     .do_norm        : normalize each row/col to go from -1 to 1
%
% OUTPUTS:
%   rho        : pairwise correlation matrix for _spikes
%   sort_idx   : indices of original _spikes inputs for sort
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho,sort_idx] = pairwiseCorr(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS
signals          =  [];
method           =  @corr;
trial_idx        =  1:length(trial_data);
cluster_order    =  false;
cluster_signals  =  false;
do_norm          =  false;
if nargin > 1, assignParams(who,params); end % overwrite defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
signals = check_signals(trial_data(1),signals);
trial_data = trial_data(trial_idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build master matrix of spiking
data = [];
for i = 1:size(signals,1)
    temp = cat(1,trial_data.(signals{i,1}));
    data = [data, temp(:,signals{i,2})];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get pairwise correlations (or other method), and replace diagonal with zeros
rho = method(data).*(-1*eye(size(data,2))+ones(size(data,2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reorder to highlight structure
if cluster_order
    if cluster_signals
        temp = cell(1,size(signals,1));
        start_idx = [0 cellfun(@length,neurons)];
        
        for i = 1:size(signals,1)
            idx = start_idx(i)+(1:neurons{i}(end));
            temp{i} = start_idx(i)+cluster_rho(rho(idx,idx));
        end
        sort_idx = cat(1,temp{:});
    else
        sort_idx = cluster_rho(rho);
    end
else
    sort_idx = 1:size(rho,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize each row/column
if do_norm
    % FIX THIS
    [temp1,temp2] = deal(zeros(size(rho)));
    temp1(rho > 0) = rho(rho > 0);
    temp1 = temp1./repmat(max(rho),size(rho,1),1);
    
    temp2(rho < 0) = rho(rho < 0);
    temp2 = -temp2./repmat(min(rho),size(rho,1),1);
    rho = temp1 + temp2;
end

% THIS IS A HACK FOR NOW REVISIT LATER
if nargout < 2
    rho = rho(sort_idx,sort_idx);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sort_idx = cluster_rho(rho)

d = squareform(pdist(rho));
c = cluster(linkage(pdist(rho),'average'),'cutoff',1);
u = unique(c);

sort_idx = [];
for i = 1:length(u)
    idx = find(c==u(i));
    
    % SORT THEM ACCORDING TO INCREASING DISTANCE
    [~,I] = sort(mean(d(idx,idx)));
    
    sort_idx = [sort_idx; idx(I)];
end

end
