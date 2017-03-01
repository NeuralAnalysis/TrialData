%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Implements Machens method for assessing dimensionality
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .signal       : which signal to work on
%       .do_smoothing : flag to smooth data
%       .kernel_SD    : SD for smoothing kernel
%       .condition    : (string) which condition to average over
%                           Default is 'target_direction'
%       .num_iter     : number of iterations (default 1000)
%       .alpha        ; what fraction of non-noise variance do you want
%                           Default is 0.95
%
% OUTPUTS:
%   dims                : how many dimensions capture alpha of non-noise variance
%   noise_eigen_prctile : 99th percentile values for the noise
%                           Note: same number of PCs as columns of signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dims,noise_eigen_prctile] = estimateDimensionality(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
signal        =  [];
sqrt_transform = false;
do_smoothing  =  false;
kernel_SD     =  0.05;
condition     =  'target_direction';
num_iter      =  1000;
alpha         =  0.95; % what fraction of non-noise variance
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(signal), error('Must provide desired signal'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbr_chs = size(trial_data(1).(signal),2);

% get number of trials to each condition
nbr_trials = cellfun(@(x) length(getTDidx(trial_data,condition,x)),num2cell(unique([trial_data.(condition)])));
tgt_idx = cellfun(@(x) getTDidx(trial_data,condition,x),num2cell(unique([trial_data.(condition)])),'uni',0);

% get PCA of smoothed, trial-averaged data
trial_data = smoothSignals(trial_data,struct('signals',signal,'do_smoothing',do_smoothing,'kernel_SD',kernel_SD,'sqrt_transform',sqrt_transform));
td = trialAverage(trial_data,condition);
[~,pca_info] = getPCA(td,struct('signals',signal));

noise_eigen = cell(1,num_iter);
[n1,n2] = size(trial_data(1).(signal));
for j = 1:num_iter
    trial_nbrs    	= randperm(min(nbr_trials),2);
    % fill the firing rate matrices
    [sfr_noise1,sfr_noise2] = deal(zeros(n1*length(tgt_idx),n2));
    for i = 1:length(tgt_idx)
        sfr_noise1((i-1)*n1+1:(i-1)*n1+n1,:) = trial_data(tgt_idx{i}(trial_nbrs(1))).(signal);
        sfr_noise2((i-1)*n1+1:(i-1)*n1+n1,:)   = trial_data(tgt_idx{i}(trial_nbrs(2))).(signal);
    end
        % calculate difference in firing rate between trials
    sfr_noise_diff  = (sfr_noise1 - sfr_noise2)/sqrt(2*min(nbr_trials));
    
    % Do PCA of the noise
    [~,~,noise_eigen{j}] = pca(sfr_noise_diff);
end


% get the amount of variance explained by the first n noise components
eigenv_noise        = cell2mat(noise_eigen);
scree_noise         = zeros(size(eigenv_noise,1),num_iter);
for j = 1:num_iter
    scree_noise(:,j) = cumsum(eigenv_noise(:,j))/sum(eigenv_noise(:,j));
end

% and turn it into a histogram
hist_x              = 0:0.01:1;
hist_scree_noise    = zeros(length(hist_x)-1,size(eigenv_noise,1));
for j = 1:size(eigenv_noise,1)
    hist_scree_noise(:,j) = histcounts(scree_noise(j,:),hist_x)/num_iter;
end

% find 99 % limit (as in Lalazar et al., PLoC Comp Biol, 2016)
noise_var_99        = zeros(1,nbr_chs);
for i = 1:nbr_chs
    noise_var_99(i) = hist_x(find(cumsum(hist_scree_noise(:,i))>0.99,1));
end
noise_eigen_prctile      = prctile(eigenv_noise,99,2);

% the cutoff is when the cumulative variance explained by the real
% eigenvalues is greater than what can be proven to be not noise
e = pca_info.eigen/sum(pca_info.eigen);
dims = find(cumsum(e)./(alpha*(1-cumsum(noise_eigen_prctile))) > 1,1,'first');

