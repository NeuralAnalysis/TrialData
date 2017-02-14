%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data,pca_info] = getPCA(trial_data, varargin)
%
% [w, mu, scores, trial_data, params] = getPCA(trial_data, params);
%   Computes PCA projection for neural data. If you request trial_data as
% a final output, will add scores to each trial for use later. Must pass in
% 'array' field for struct. Note that this can be a cell of strings to pool
% data from multiple arrays. In this case, the rows of w will be as if you
% concatenated the two arrays together in the order they were provided.
%
% trial_data = getPCA(trial_data, w, mu, params);
%   Uses an old w and mu from previous getPCA call to add scores to
% trial_data as ARRAY_pca. Params is still needed to specify the array or
% if you want smoothing, etc.
%
%   NOTE: always de-means! Theoretically this could be modified to take a
% parameter to skip the de-meaning.
%
% INPUTS:
%   trial_data : the struct
%   params     : struct containing parameters
%     .array          : which units (can be cell array with multiple)
%     .trial_idx      : which trials to use (default: all)
%     .neurons        : which neurons to use (default: all) Note: for multiple arrays,
%                       neurons should be cell array with indices for each array
%     .sqrt_transform : flag to square root transform spikes (default: true)
%     .do_smoothing   : flag to convolve spikes with gaussian (default: true)
%     .kernel_SD      : kernel s.d. for smoothing (default: 2*bin_size)
%     .trial_avg      : flag to trial average (requires condition input) (default: false)
%     .trial_avg_cond : (string/cell) which conditions to average over
%     .do_plot        : flag to make scree plot (default: false)
%
% OUTPUTS:
%   w          : weight matrix for PCA projections
%   scores     : scores for the PCs
%   eigen      : eigenvalues for PC ranking
%   mu         : mean for each input (for de-meaning later)
%                  NOTE: if you use w later, you MUST demean using mu!!!
%   trial_data : old struct with added field for scores for each trial
%                   NOTE: if passing in old w and mu, only returns this
%
% EXAMPLES:
%   e.g. to compute covariance matrix
%       [w,mu,~,~,params] = getPCA(trial_data, struct('array','M1','bin_size',0.01));
%   e.g. to add scores to trial_data later using the above output
%       trial_data = getPCA(trial_data, w, mu, params);
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,pca_info] = getPCA(trial_data, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if length(varargin) == 1 % only params is provided
    new_pca = true;
    params = varargin{1};
else
    new_pca = false;
    if nargout > 1, error('When using old PCA, will only output trial_data'); end
    if length(varargin) == 3 % provided cov matrix and mu
        w = varargin{1};
        mu = varargin{2};
        params = varargin{3};
    else
        error('Incorrect number of inputs');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get defaults or process parameters
if ~exist('params','var'), params = struct(); end
if isfield(params,'array'),array = params.array; else, error('Need to specify array'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
trial_idx       =  1:length(trial_data);
neurons         =  [];
sqrt_transform  =  true;
do_smoothing    =  true;
kernel_SD       =  0.05;
trial_avg       =  false;
trial_avg_cond  =  {};
do_plot         =  false;
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if trial_avg && isempty(trial_avg_cond), error('Must provide conditions to average trials over.'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process and prepare inputs
if ~iscell(array), array = {array}; end
if isempty(neurons)
    for i = 1:length(array), neurons{i} = 1:size(trial_data(1).([array{i} '_spikes']),2); end
end
if ~iscell(neurons), neurons = {neurons}; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop along trials to square root transform and smooth if desired
td = smoothSpikes(trial_data(trial_idx),params);
if trial_avg, td = trialAverage(td,trial_avg_cond); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concatenate specified trials
fr = [];
for i = 1:length(array)
    temp_fr = cat(1,td.([array{i} '_spikes']));
    fr = [fr, temp_fr(:,neurons{i})];
end
% get the time points that separate each trial later
trial_markers = [1,1+cumsum(cellfun(@(x) size(x,1),{td.pos}))];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build PCA model for M1
if new_pca
    % compute PCA
    [w, scores, eigen,~,~,mu] = pca(fr,'Algorithm','svd');
    
    if do_plot
        figure,
        subplot(2,1,1);
        bar(eigen/sum(eigen));
        axis('tight');
        xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
        set(gca,'TickDir','out'),set(gca,'FontSize',14);
        xlim([0 size(fr,2)+1])
        
        subplot(2,1,2);
        plot(cumsum(eigen/sum(eigen)),'linewidth',3,'marker','d'),
        xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
        set(gca,'TickDir','out'),set(gca,'FontSize',14);
        xlim([0 size(fr,2)+1])
        ylim([0 1])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add scores to trial_data
if trial_avg, trial_data = td; end % THIS IS A HACK FOR NOW
for trial = 1:length(trial_data)
    idx = trial_markers(trial):trial_markers(trial+1)-1;
    trial_data(trial).([[array{:}] '_pca']) = (fr(idx,:)-repmat(mu,length(idx),1))*w;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package up outputs
if new_pca && nargout == 2
    pca_params = struct('array',array,'sqrt_transform',sqrt_transform,'do_smoothing',do_smoothing,'bin_size',bin_size,'kernel_SD',kernel_SD);
    pca_info = struct('w',w,'mu',mu,'scores',scores,'eigen',eigen,'params',pca_params);
end