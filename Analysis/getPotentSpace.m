%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data, pca_info] = getPotentSpace(trial_data, params)
%
%   This will compute the geometry of the potent and null spaces for two
% sets of neurons. Calls getPCA for neural data processing.
%
% TO DO:
%   - use arbitrary signals like EMG or velocity as output space
%
% INPUTS:
%   trial_data : the struct
%   params    : struct containing parameters
%     .in_array    : (string) name of input array
%     .out_array   : (string) name of output array
%     .in_dims     : dimensionality of input space
%     .out_dims    : dimensionality of output space
%     .in_neurons  : list of neurons to include for in array (default: all)
%     .out_neurons : list of neurons to include for out array (default: all)
%     .use_trials  : vector list of trial indices to use for spaces
%         NOTE: if adds scores to trial_data, will add scores for all, not only use_trials
%     (Also include any parameters desired for getPCA, as params gets
%     passed to that function call as well)
%
% OUTPUTS:
%   trial_data : the struct with potent/null fields added
%                NOTE: processed by getPCA too (e.g. trial averaging, etc)
%   pca_results : struct with PCA output
%       .V_potent : potent space basis vectors
%       .V_null   : null space basis vectors
%       .w_in     : weights for PC space of input array
%       .w_out    : weights for PC space of output array
%       .mu_in    : means of input array
%       .mu_out   : means of output array
%                NOTE: as with PCA, the mu are essential for V/w later
%       .params   : parameter struct
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [td,pca_info] = getPotentSpace(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
in_array     =  [];
out_array    =  [];
in_dims      =  [];
out_dims     =  [];
in_neurons   =  [];
out_neurons  =  [];
use_trials   =  1:length(trial_data);
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(in_array), error('Need to specify input array name'); end
if isempty(out_array), error('Need to specify output array name'); end
if isempty(in_dims), error('Need to specify input dimensionality'); end
if isempty(out_dims), error('Need to specify output dimensionality'); end
if isempty(in_neurons), in_neurons = 1:size(trial_data(1).([in_array '_spikes']),2); end
if isempty(out_neurons), out_neurons = 1:size(trial_data(1).([out_array '_spikes']),2); end

pca_params = params;

% get output PC space
pca_params.arrays = out_array;
pca_params.neurons = out_neurons;
[~,pca_info] = getPCA(trial_data(use_trials),pca_params);
w_out = pca_info.w;
mu_out = pca_info.mu;
score_out = pca_info.scores;

% get input PC space
pca_params.arrays = in_array;
pca_params.neurons = in_neurons;
[td,pca_info] = getPCA(trial_data(use_trials),pca_params);
w_in = pca_info.w;
mu_in = pca_info.mu;
score_in = pca_info.scores;
pca_params = pca_info.params;

y = score_out(:,1:out_dims);
x = score_in(:,1:in_dims);
% find the model
W = zeros( size(y,2), size(x,2) );
for i = 1:size(y,2)
    %[b_pc, ~, ~, ~, stats_this] = regress(y(:,i),x);
    b_pc = x\y(:,i);
    % fill MIMO matrix W
    W(i,:) = b_pc';
end

% do SVD of weights to get potent/null spaces
[U, S, V]                   = svd( W );
% The output potent spaces is defined by the first m columns of V', where m
% is the number of dimensions of the output
V_potent                    = V(1:size(y,2),:)';
V_null                      = V(size(y,2)+1:end,:)';

% package up PCA weights, etc.
pca_info = struct(        ...
    'V_potent', V_potent, ...
    'V_null',   V_null,   ...
    'w_in',     w_in,     ...
    'w_out',    w_out,    ...
    'mu_in',    mu_in,    ...
    'mu_out',   mu_out,   ...
    'params',   pca_params);

for trial = 1:length(td)
    temp = td(trial).([in_array '_pca']);
    temp = temp(:,1:in_dims);
    td(trial).potent = temp*V_potent;
    td(trial).null = temp*V_null;
end



