%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho,sort_idx] = pairwiseCorr(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Compute pairwise correlations between neurons. Can order them to show
% structure, a la Cunningham.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFAULTS
arrays          =  [];
trial_idx       =  1:length(trial_data);
neurons         =  [];
cluster_order   =  false;
array_cluster   =  true;
if nargin > 1, assignParams(who,params); end % overwrite defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if ~iscell(arrays), arrays = {arrays}; end
if isempty(neurons)
    for i = 1:length(arrays), neurons{i} = size(trial_data(1).([arrays{i} '_spikes']),2); end
end
if ~iscell(neurons), neurons = {neurons}; end
if isempty(arrays), error('Must give array'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build master matrix of spiking
fr = [];
for array = 1:length(arrays)
    temp = cat(1,td.([arrays{array} '_spikes']));
    fr = [fr, temp(:,neurons{array})];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get pairwise correlations, and replace diagonal with zeros
rho = corr(fr).*(-1*eye(size(fr,2))+ones(size(fr,2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reorder to highlight structure

% if cluster_order
%     sort_idx = cell(1,length(arrays));
%     if array_cluster
%         for array = 1:length(arrays)
%             sort_idx{array} = neurons{array}(idx);
%         end
%     else
%         for array = 1:length(arrays)
%             
%             sort_idx{array} = neurons{array}(idx);
%         end
%     end
%     rho = rho(idx,idx);
% else
%     sort_idx = neurons;
% end

