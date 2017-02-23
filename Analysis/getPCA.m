%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data,pca_info] = getPCA(trial_data, varargin)
%
% [trial_data,pca_info] = getPCA(trial_data, params);
%   Computes PCA projection for time-varying data. Will add scores to each trial
% for use later. Must pass in 'signal' field for struct. Note that this can
% be a cell of strings to pool data from multiple arrays. In this case, the
% rows of w will be as if you concatenated the columns of signals together in the
% order they were provided.
%
% trial_data = getPCA(trial_data, w, params);
%   Uses an old w from previous getPCA call to add scores to
% trial_data as SIGNAL_pca. Params is still needed to specify the array or
% if you want smoothing, etc.
%
%   NOTE: centers data by default! Thus to reconstruct scores you need the
%   means of each signal
%
% INPUTS:
%   trial_data : the struct
%   params     : struct containing parameters
%     .signals        : which field names (can be cell with multiple)
%                           NOTE: Defaults to using all the _spikes fields
%     .trial_idx      : which trials to use (default: all)
%     .signal_idx     : which columns to use (default: all) Note: for multiple signals,
%                       neurons should be cell array with indices for each signal
%     .sqrt_transform : flag to square root transform spikes (default: false)
%     .do_smoothing   : flag to convolve spikes with gaussian (default: false)
%     .kernel_SD      : kernel s.d. in s for smoothing (default: 0.05)
%     .trial_avg_cond : (string/cell) which conditions to average over
%                          NOTE: if empty or not passed in, won't trial average
%     .do_plot        : flag to make scree plot (default: false)
%
% OUTPUTS:
%   trial_data : old struct with added field for scores for each trial
%                   NOTE: if passing in old w, only returns this
%   pca_info   : struct of PCA information
%     .w          : weight matrix for PCA projections
%     .scores     : scores for the PCs
%     .eigen      : eigenvalues for PC ranking
%     .params     : the parameters used for this analysis
%
% EXAMPLES:
%   e.g. to compute covariance matrix
%       [~,pca_info] = getPCA(trial_data, struct('signals','M1_spikes'));
%       w = pca_info.w;
%   e.g. to add scores to trial_data later using the above output
%       trial_data = getPCA(trial_data, w, params);
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
elseif  isempty(varargin)
    new_pca = true;
    params = struct();
else
    new_pca = false;
    if nargout > 1, error('When using old PCA, will only output trial_data'); end
    if length(varargin) == 3 % provided cov matrix
        w = varargin{1};
        params = varargin{2};
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
signals         =  getTDfields(trial_data,'spikes');
signal_idx      =  [];
sqrt_transform  =  false;
do_smoothing    =  false;
kernel_SD       =  0.05;
trial_avg_cond  =  {};
do_plot         =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters you can change that aren't described in header
pca_algorithm   = 'svd'; % algorithm for PCA
pca_centered    = true;  % whether to center data
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process and prepare inputs
if isempty(trial_avg_cond), trial_avg = false; else, trial_avg = true; end
if ~iscell(signals), signals = {signals}; end
if isempty(signal_idx)
    signal_idx = cell(1,length(signals));
    for i = 1:length(signals), signal_idx{i} = 1:size(trial_data(1).(signals{i}),2); end
end
if ~iscell(signal_idx), signal_idx = {signal_idx}; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% square root transform and smooth if desired
smooth_params = struct( ...
    'signals',        signals, ...
    'sqrt_transform', sqrt_transform, ...
    'do_smoothing',   do_smoothing, ...
    'kernel_SD',      kernel_SD);
td = smoothSignals(trial_data(trial_idx),smooth_params);
if trial_avg, td = trialAverage(td,trial_avg_cond); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concatenate specified trials
data = [];
for i = 1:length(signals)
    temp_data = cat(1,td.(signals{i}));
    data = [data, temp_data(:,signal_idx{i})];
end
% get the time points that separate each trial later
fn_time = getTDfields(td,'time');
trial_markers = [1,1+cumsum(cellfun(@(x) size(x,1),{td.(fn_time{1})}))];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build PCA model for M1
if new_pca
    % compute PCA
   [w, scores, eigen] = pca(data,'Algorithm',pca_algorithm,'Centered',pca_centered);
    
    if do_plot
        figure,
        subplot(2,1,1);
        bar(eigen/sum(eigen));
        axis('tight');
        xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
        set(gca,'Box','off','TickDir','out','FontSize',14);
        xlim([0 size(data,2)+1])
        
        subplot(2,1,2);
        plot(cumsum(eigen/sum(eigen)),'linewidth',3,'marker','d'),
        xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
        set(gca,'Box','off','TickDir','out','FontSize',14);
        xlim([0 size(data,2)+1])
        ylim([0 1])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package up outputs
if new_pca
    pca_params = struct( ...
        'signals',signals, ...
        'trial_idx',trial_idx, ...
        'signal_idx',{signal_idx}, ...
        'sqrt_transform',sqrt_transform, ...
        'do_smoothing',do_smoothing, ...
        'kernel_SD',kernel_SD, ...
        'trial_avg_cond',trial_avg_cond);
    pca_info = struct('w',w,'scores',scores,'eigen',eigen,'params',pca_params);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add scores to trial_data
if pca_centered, mu = mean(data,1); else mu = zeros(1,size(data,2)); end
signals = cellfun(@(x) strrep(x,'_spikes',''),signals,'uni',0);
for trial = 1:length(trial_data)
    idx = trial_markers(trial):trial_markers(trial+1)-1;
    trial_data(trial).([[signals{:}] '_pca']) = (data(idx,:)-repmat(mu,length(idx),1))*w;
end

