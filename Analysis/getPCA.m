%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data,pca_info] = getPCA(trial_data, varargin)
%
% [trial_data,pca_info] = getPCA(trial_data, params);
%   Computes PCA projection for neural data. Will add scores to each trial
% for use later. Must pass in 'array' field for struct. Note that this can
% be a cell of strings to pool data from multiple arrays. In this case, the
% rows of w will be as if you concatenated the two arrays together in the
% order they were provided.
%
% TO DO:
%   - make more generalizable (any input e.g. EMG instead of arrays only)
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
%     .arrays         : which arrays (can be cell with multiple)
%                           NOTE: Defaults to using all the _spikes fields
%     .trial_idx      : which trials to use (default: all)
%     .neurons        : which neurons to use (default: all) Note: for multiple arrays,
%                       neurons should be cell array with indices for each array
%     .sqrt_transform : flag to square root transform spikes (default: false)
%     .do_smoothing   : flag to convolve spikes with gaussian (default: false)
%     .kernel_SD      : kernel s.d. in s for smoothing (default: 0.05)
%     .trial_avg_cond : (string/cell) which conditions to average over
%                          NOTE: if empty or not passed in, won't trial average
%     .do_plot        : flag to make scree plot (default: false)
%
% OUTPUTS:
%   trial_data : old struct with added field for scores for each trial
%                   NOTE: if passing in old w and mu, only returns this
%   pca_info   : struct of PCA information
%     .w          : weight matrix for PCA projections
%     .scores     : scores for the PCs
%     .eigen      : eigenvalues for PC ranking
%     .mu         : mean for each input (for de-meaning later)
%                     NOTE: if you use w later, you MUST demean using mu!!!
%     .params     : the parameters used for this analysis
%
% EXAMPLES:
%   e.g. to compute covariance matrix
%       [~,pca_info] = getPCA(trial_data, struct('arrays','M1'));
%       w = pca_info.w;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
trial_idx       =  1:length(trial_data);
arrays          =  getTDfields(trial_data,'arrays');
neurons         =  [];
sqrt_transform  =  false;
do_smoothing    =  false;
kernel_SD       =  0.05;
trial_avg_cond  =  {};
do_plot         =  false;
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process and prepare inputs
if isempty(trial_avg_cond), trial_avg = false; else, trial_avg = true; end
if ~iscell(arrays), arrays = {arrays}; end
if isempty(neurons)
    for i = 1:length(arrays), neurons{i} = 1:size(trial_data(1).([arrays{i} '_spikes']),2); end
end
if ~iscell(neurons), neurons = {neurons}; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% square root transform and smooth if desired
smooth_params = struct( ...
    'sqrt_transform',sqrt_transform, ...
    'do_smoothing',do_smoothing, ...
    'kernel_SD',kernel_SD);
td = smoothSpikes(trial_data(trial_idx),smooth_params);
if trial_avg, td = trialAverage(td,trial_avg_cond); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concatenate specified trials
fr = [];
for i = 1:length(arrays)
    temp_fr = cat(1,td.([arrays{i} '_spikes']));
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
    trial_data(trial).([[arrays{:}] '_pca']) = (fr(idx,:)-repmat(mu,length(idx),1))*w;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package up outputs
if new_pca && nargout == 2
    pca_params = struct( ...
        'arrays',arrays, ...
        'trial_idx',trial_idx, ...
        'neurons',{neurons}, ...
        'sqrt_transform',sqrt_transform, ...
        'do_smoothing',do_smoothing, ...
        'kernel_SD',kernel_SD, ...
        'trial_avg_cond',trial_avg_cond);
    pca_info = struct('w',w,'mu',mu,'scores',scores,'eigen',eigen,'params',pca_params);
end