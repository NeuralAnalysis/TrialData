function [V_potent, V_null, w_cov, w_pred, mu_cov, mu_pred] = getPotentSpace(trial_data,cov_array,pred_array,params)

bin_size = params.bin_size;
pca_dims = params.pca_dims;

num_pred_neurons = size(trial_data(1).([pred_array '_spikes']),2);
num_cov_neurons = size(trial_data(1).([cov_array '_spikes']),2);

if 0
    disp('USING BASELINE FOR POTENT');
    [w_cov,score_cov, ~, mu_cov] = getPCA(trial_data,struct('array',cov_array,'bin_size',params.dt,'trial_idx',find(getTDidx(trial_data,'epoch','BL')),'neurons',1:num_cov_neurons));
else
    [w_cov,score_cov, ~, mu_cov] = getPCA(trial_data,struct('array',cov_array,'bin_size',params.dt,'trial_idx',1:length(trial_data),'neurons',1:num_cov_neurons));
end
params.pca_w = w_cov;
% do SVD to get null/potent spaces

% get M1 PC space
if 0
    [w_pred, score_pred, ~, mu_pred] = getPCA(trial_data,struct('array',pred_array,'bin_size',params.dt,'trial_idx',find(getTDidx(trial_data,'epoch','BL')),'neurons',1:num_pred_neurons));
else
    [w_pred, score_pred, ~, mu_pred] = getPCA(trial_data,struct('array',pred_array,'bin_size',params.dt,'trial_idx',1:length(trial_data),'neurons',1:num_pred_neurons));
end
y = score_pred(:,pca_dims.(pred_array));
x = score_cov(:,pca_dims.(cov_array));
% find the model
W = zeros( size(y,2), size(x,2) );
for i = 1:size(y,2)
    %[b_pc, ~, ~, ~, stats_this] = regress(y(:,i),x);
    b_pc = x\y(:,i);
    % fill MIMO matrix W
    W(i,:) = b_pc';
end
% do SVD of weights
[U, S, V]                   = svd( W );
% The output potent spaces is defined by the first m columns of V', where m
% is the number of dimensions of the output
V_potent                    = V(1:size(y,2),:)';
V_null                      = V(size(y,2)+1:end,:)';
