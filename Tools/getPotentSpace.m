%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [V_potent, V_null, w_in, w_out, mu_cov, mu_pred] = 
%               getPotentSpace(trial_data,cov_array,pred_array,params)
%
%   This will compute the geometry of the potent and null spaces for two
% sets of neurons. This may be a bit specialized for Matt's needs at the
% moment. Calls getPCA.
%
% TO DO:
%   - add scores from each space to trial_data
%   - give option to select trial subsets (and get rid of do_bl thing)
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
%
% OUTPUTS:
%   V_potent : projection matrix for potent space
%   V_null   : projection matrix for null space
%   w_in     : weights for PC space of input array
%   w_out    : weights for PC space of output array
%   mu_in    : means of input array
%   mu_out   : means of output array
%                NOTE: as with PCA, the mu are essential for V/w later
% 
% Written by Matt Perich. Updated Feb 2017.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V_potent, V_null, w_in, w_out, mu_in, mu_out] = getPotentSpace(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(params,'in_array'), in_array = params.in_array; else, error('Need to specify input array name'); end
if isfield(params,'out_array'), out_array = params.out_array; else, error('Need to specify output array name'); end
if isfield(params,'in_dims'), in_dims = params.in_dims; else, error('Need to specify input dimensionality'); end
if isfield(params,'out_dims'), out_dims = params.out_dims; else, error('Need to specify output dimensionality'); end
if isfield(params,'in_neurons'), in_neurons = params.in_neurons; else, in_neurons = 1:size(trial_data(1).([in_array '_spikes']),2); end
if isfield(params,'out_neurons'), out_neurons = params.out_neurons; else, out_neurons = 1:size(trial_data(1).([out_array '_spikes']),2); end

% dumb hard coded thing for matt's personal use. SHOULD BE REMOVED
do_bl = false;

if do_bl
    disp('USING BASELINE FOR POTENT');
    [w_in,score_in, ~, mu_in] = getPCA(trial_data,struct('array',in_array,'bin_size',params.dt,'trial_idx',find(getTDidx(trial_data,'epoch','BL')),'neurons',in_neurons));
else
    [w_in,score_in, ~, mu_in] = getPCA(trial_data,struct('array',in_array,'bin_size',params.dt,'trial_idx',1:length(trial_data),'neurons',in_neurons));
end
params.pca_w = w_in;

% get output PC space
if do_bl
    [w_out, score_out, ~, mu_out] = getPCA(trial_data,struct('array',out_array,'bin_size',params.dt,'trial_idx',find(getTDidx(trial_data,'epoch','BL')),'neurons',out_neurons));
else
    [w_out, score_out, ~, mu_out] = getPCA(trial_data,struct('array',out_array,'bin_size',params.dt,'trial_idx',1:length(trial_data),'neurons',out_neurons));
end

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
