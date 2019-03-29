%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Implements Machens method for assessing dimensionality.
%
% IMPORTANT: Must be done on single-trial data. No trialAverage!
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .signals      : which signal(s) to work on
%       .condition    : (string) which condition to average over
%                           Default is 'target_direction'
%       .use_trials   : which trials to use (default is all)
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
signals         =  [];
condition       =  'target_direction';
use_trials      =  1:length(trial_data);
num_iter        =  1000;
alpha           =  0.95; % what fraction of non-noise variance
trim_idx        =  {};   % can trim data in here {'idx',val;'idx',val}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some undocumented parameters that you can overwrite if you need to
pca_algorithm   =  'svd'; % algorithm for PCA
pca_economy     =  false; % if num samples < degrees of freedom, will pad with zeros to keep output same size as degrees of freedom
pca_centered    =  true;  % whether to center data
verbose         =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data  =  check_td_quality(trial_data);

if isempty(signals), error('Must provide desired signal'); end
signals = check_signals(trial_data(1),signals);
if iscell(use_trials) % likely to be meta info
    use_trials = getTDidx(trial_data,use_trials{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trial_data = trial_data(use_trials);

% Process the desired signals
if size(signals,1) > 1
    % this is a hack, but make a new field that is the concatenated signals
    signal = [signals{:,1}];
    for trial = 1:length(trial_data)
        trial_data(trial).(signal) = [];
        for i = 1:size(signals,1)
            temp = trial_data(trial).(signals{i,1});
            trial_data(trial).(signal) = [trial_data(trial).(signal), temp(:,signals{i,2})];
        end
    end
else
    signal = signals{1,1};
    for trial = 1:length(trial_data)
        temp = trial_data(trial).(signal);
        trial_data(trial).(signal) = temp(:,signals{1,2});
    end
end

nbr_chs = size(trial_data(1).(signal),2);

% get number of trials to each condition
nbr_trials = cellfun(@(x) length(getTDidx(trial_data,condition,x)),num2cell(unique([trial_data.(condition)])));
tgt_idx = cellfun(@(x) getTDidx(trial_data,condition,x),num2cell(unique([trial_data.(condition)])),'uni',0);

if any(nbr_trials < 2)
    error('Too few trials for one or more conditions.');
end

if ~isempty(trim_idx)
    trial_data = trimTD(trial_data,trim_idx(1,:),trim_idx(2,:));
end


% get PCA of smoothed, trial-averaged data
td = trialAverage(trial_data,condition);

if size(cat(1,td.(signal)),1) < size(cat(1,td.(signal)),2)
    warning('Number of total datapoints across trials is smaller than the total degrees of freedom! Be careful...');
end

[~,pca_info] = dimReduce(td,struct('signals',signal,'algorithm','pca','pca_economy',pca_economy,'pca_algorithm',pca_algorithm,'pca_centered',pca_centered));

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
    
    % This will makea  plot of a single neuron and the noise etc
    % % %     [sfr_noise1,sfr_noise2] = deal(zeros(n1,n2));
    % % %     for i = 1%:length(tgt_idx)
    % % %         sfr_noise1((i-1)*n1+1:(i-1)*n1+n1,:) = trial_data(tgt_idx{i}(trial_nbrs(1))).(signal);
    % % %         sfr_noise2((i-1)*n1+1:(i-1)*n1+n1,:)   = trial_data(tgt_idx{i}(trial_nbrs(2))).(signal);
    % % %     end
    % % %         % calculate difference in firing rate between trials
    % % %     sfr_noise_diff  = (sfr_noise1 - sfr_noise2)/sqrt(2*min(nbr_trials));
    % % %     figure;
    % % %     subplot(121); hold all;
    % % %     for i = 1:nbr_trials, temp = trial_data(tgt_idx{1}(i)).(signal); plot(temp(:,28),'k'); end;
    % % %     temp = td(1).(signal);
    % % %     plot(temp(:,28),'r-','LineWidth',2);
    % % %     set(gca,'Box','off','TickDir','out','FontSize',14);
    % % %     subplot(122); hold all;
    % % %     for i = 1:size(sfr_noise_diff,2)
    % % %         if i == 28
    % % %             plot(sfr_noise_diff(:,i),'r','LineWidth',2);
    % % %         else
    % % %             plot(sfr_noise_diff(:,i),'k');
    % % %         end
    % % %     end
    % % %     set(gca,'Box','off','TickDir','out','FontSize',14);
    
    
    
    % Do PCA of the noise
    [~,~,noise_eigen{j}] = pca(sfr_noise_diff,'Algorithm',pca_algorithm,'Centered',pca_centered,'Economy',pca_economy);
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

% % % ANOTHER PLOT FOR ILLUSTRATING THE PROCESS
% % figure; hold all;
% % plot(hist_scree_noise);
% % set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 100],'YLim',[0 1]);
% % xlabel('% Variance Explained');
% % ylabel('Density');

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


% %% Make some more fun plots to show the method
% % % plot noise cumsum and threshold
% % figure;
% % subplot(211); hold all;
% % plot(cumsum(noise_eigen_prctile),'LineWidth',2);
% % % plot([1,length(noise_eigen_prctile)],[alpha,alpha],'k--','LineWidth',1);
% % plot([1,length(noise_eigen_prctile)],[(alpha*(sum(noise_eigen_prctile))),(alpha*(sum(noise_eigen_prctile)))],'k--','LineWidth',1);
% % axis('tight');
% % set(gca,'Box','off','TickDir','out','FontSize',14);
% % xlabel('Nbr Noise Dims');
% % ylabel('Cum Var Explained');
% %
% % subplot(212); hold all;
% % bar(e);
% % plot([1,length(noise_eigen_prctile)],[(alpha*(sum(noise_eigen_prctile))),(alpha*(sum(noise_eigen_prctile)))],'k--','LineWidth',1);
% % set(gca,'Box','off','TickDir','out','FontSize',14);
% % axis('tight');
% %
% % xlabel('Nbr Neural Dims');
% % ylabel('Var Explained');
% % V = axis;
% % plot([dims,dims],V(3:4),'r-','LineWidth',2);

