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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dims,noise_eigen_99] = estimateDimensionality(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
signal        =  [];
sqrt_transform = true;
do_smoothing  =  true;
kernel_SD     =  0.05;
condition     =  'target_direction';
nbr_iter      =  1000;
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
trial_data = smoothSignal(trial_data,struct('signals',signal,'do_smoothing',do_smoothing,'kernel_SD',kernel_SD,'sqrt_transform',sqrt_transform));
td = trialAverage(trial_data,condition);
[td,pca_info] = getPCA(td,struct('signals',signal));

for j = 1:nbr_iter
    trial_nbrs    	= randperm(min(nbr_trials),2);
    % fill the firing rate matrices
    sfr_noise    	= [];
    sfr_noise2     	= [];
    for i = 1:length(tgt_idx)
        sfr_noise   = [sfr_noise; trial_data(tgt_idx{i}(trial_nbrs(1))).(signal)];
        sfr_noise2  = [sfr_noise2; trial_data(tgt_idx{i}(trial_nbrs(2))).(signal)];
    end
        % calculate difference in firing rate between trials
    sfr_noise_diff  = (sfr_noise - sfr_noise2)/sqrt(2*min(nbr_trials));
    
    % Do PCA of the noise
    [~,~,noise_eigen{j},~,~,~] = pca(sfr_noise_diff);
end


% get the amount of variance explained by the first n noise components
eigenv_noise        = cell2mat(noise_eigen);
scree_noise         = zeros(size(eigenv_noise,1),nbr_iter);
for j = 1:nbr_iter
    scree_noise(:,j) = cumsum(eigenv_noise(:,j))/sum(eigenv_noise(:,j));
end

% and turn it into a histogram
hist_x              = 0:0.01:1;
hist_scree_noise    = zeros(length(hist_x)-1,size(eigenv_noise,1));
for j = 1:size(eigenv_noise,1)
    hist_scree_noise(:,j) = histcounts(scree_noise(j,:),hist_x)/nbr_iter;
end

% find mean noise variances
mean_noise_var      = mean(cumsum(eigenv_noise),2)/mean(sum(eigenv_noise),2);
mean_noise_eigen     = mean(eigenv_noise,2);

% find 99 % limit (as in Lalazar et al., PLoC Comp Biol, 2016)
noise_var_99        = zeros(1,nbr_chs);
for i = 1:nbr_chs
    noise_var_99(i) = hist_x(find(cumsum(hist_scree_noise(:,i))>0.99,1));
end
noise_eigen_99      = prctile(eigenv_noise,99,2);


e = pca_info.eigen/sum(pca_info.eigen);
dims = find(cumsum(e)./(alpha*(1-cumsum(noise_eigen_99))) > 1,1,'first');

