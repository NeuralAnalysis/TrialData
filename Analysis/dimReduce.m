%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data,info_out] = dimReduce(trial_data, params)
%
% [trial_data,info_out] = dimReduce(trial_data, params);
%   Dimensionality-reduction for time-varying data. Will add scores to each trial
% for use later. Must pass in 'signals' field for struct. Note that this can
% be a cell of strings to pool data from multiple arrays. In this case, the
% rows of w will be as if you concatenated the columns of signals together in the
% order they were provided.
%
% trial_data = dimReduce(trial_data, info_out);
%   Uses an old w from previous dimReduce call to add scores to
% trial_data as SIGNALS{:}_ALGORITHM. (e.g. for M1 signals processed with
% PCA, the fieldname is 'M1_pca').
%
% Several algorithms are supported. Check the paramter list, since each can
% have specific options:
%   'pca'   : vanilla PCA, based on Matlab's built-in pca function
%   'ppca   : probabalistc PCA, based on Matlab's built-in ppca function
%   'fa'    : factor analysis, based on Matlab's factoran
%       NOTE: for PPCA/FA, must provide num_dims as input
%
% NOTE: centers data by default! Thus to reconstruct scores you need the
%       means of each signal.
%
% A hint: if you want the default parameters, you can pass in just the
% signals that you want instead of a params struct.
%   e.g., trial_data = dimReduce(trial_data,'M1_spikes');
%
% INPUTS:
%   trial_data : the struct
%   params     : struct containing parameters
%     .signals        : (cell) which signals. Two options:
%                           1) {'NAME1','NAME2',etc}
%                           2) {'NAME1',idx; 'NAME2',idx; etc}
%                                   Here idx is which columns to use
%                                   Note: can use 'all' as idx for all
%     .trial_idx      : which trials to use for computing the low-D space (default: all)
%                         Note: when adding projections, does so for all
%                         trials in trial_data, not just trial_idx
%     .num_dims       : how many dimensions (e.g. for PPCA, FA). For PCA,
%                        do nothing and it returns same dimensionality as
%                        input, or specify a value if you like.
%     .do_plot        : flag to make scree plot (default: false)
%
% OUTPUTS:
%   trial_data : old struct with added field for scores for each trial
%                   NOTE: if passing in old w, only returns this
%   info_out   : struct of information. Contents can vary based on algorithm
%     .w          : weight matrix for projections
%     .scores     : scores for the components
%     .eigen      : eigenvalues for PC ranking
%     .params     : the parameters used for this analysis
%
% EXAMPLES:
%   e.g. to compute covariance matrix
%       [~,info_out] = dimReduce(trial_data, struct('signals','M1_spikes'));
%   e.g. to add scores to trial_data later using the above output
%       trial_data = dimReduce(trial_data, params);
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,info_out] = dimReduce(trial_data, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
algorithm       = 'pca';
signals         =  getTDfields(trial_data,'spikes');
use_trials      =  1:length(trial_data);
num_dims        =  []; % how many dimensions, needed for PPCA, FA, etc
do_plot         =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters you can change that aren't described in header
sig_name         = '';    % output will be in field "SIG_NAMES_ALGORITHM". Defaults to concatenated names of signals
sqrt_transform   = false; % square root transform before reduction (projections don't have it) 
do_smoothing     = false; % will smooth before dim reduction  (trial_data projections are unsmoothed)
kernel_SD        = 0.05;  %   gaussian kernel s.d. for smoothing
pca_algorithm    = 'svd'; % algorithm for PCA
pca_economy      = false; % if num samples < degrees of freedom, will pad with zeros to keep output same size as degrees of freedom
pca_centered     = true;  % whether to center data
fa_orthogonalize = true; % whether to orthogonalize the projections
fa_rotate        = 'varimax';
add_proj_to_td   = true;  % whether to add projections
recenter_for_proj = false; % whether to recenter data before projecting into PC space
w                 = [];    % w is used to know if params was info_out (e.g. whether to recompute space)
mu                = [];    % mu is the mean from fitting, only filled if info_out is passed in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if ~isstruct(params)
        % assume the signals were provided
        signals = params;
    else
        assignParams(who,params); % overwrite parameters
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process and prepare inputs
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end
signals = check_signals(trial_data(1),signals);
if iscell(use_trials) % likely to be meta info
    use_trials = getTDidx(trial_data,use_trials{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if info_out not already sent in...build model
if isempty(w)
    td = trial_data(use_trials);
    if sqrt_transform
        td = sqrtTransform(td,signals);
    end
    if do_smoothing
        td = smoothSignals(td,struct('signals',{signals},'kernel_SD',kernel_SD));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % concatenate specified trials
    data = [];
    for i = 1:size(signals,1)
        temp_data = cat(1,td.(signals{i,1}));
        data = [data, temp_data(:,signals{i,2})];
    end
    clear td;
    
    if size(data,1) < size(data,2)
        warning('Number of total datapoints across trials is smaller than the total degrees of freedom! Be careful...');
    end
    
    % compute reduced space
    switch lower(algorithm)
        case 'pca'
            [w, scores, eigen,~,~,mu] = pca(data,'Algorithm',pca_algorithm,'Centered',pca_centered,'Economy',pca_economy);
            stats = [];
        case 'ppca'
            if isempty(num_dims), error('Must specify number of dimensions when using PPCA'); end
            [w, scores, eigen, mu] = ppca(data,num_dims);
            stats = [];
        case 'fa'
            if isempty(num_dims), error('Must specify number of dimensions when using FA'); end
            [w, eigen, ~, stats, scores] = factoran(data,num_dims,'Rotate',fa_rotate);
            mu = zeros(1,size(data,2));
    end
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Package up outputs
    params = struct( ...
        'signals',{signals}, ...
        'trial_idx',use_trials);
    info_out = struct('w',w,'scores',scores,'eigen',eigen,'mu',mu,'signals',{signals},'params',params,'sig_name',sig_name,'stats',stats);
else
    info_out = params;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add scores to trial_data
if add_proj_to_td
    if isempty(sig_name)
        sig_name = cellfun(@(x) strrep(x,'_spikes',''),signals(:,1),'uni',0);
    else
        if ~iscell(sig_name), sig_name = {sig_name}; end
    end
    n_signals = cellfun(@(x) length(x),signals(:,2));
    
    if recenter_for_proj
        if pca_centered
            if ~strcmpi(algorithm,'fa')
                mu = mean(get_vars(trial_data,signals),1);
            else
                warning('Need to look into how to handle centering with FA');
                mu = zeros(1,sum(n_signals));
            end
        else
            mu = zeros(1,sum(n_signals));
        end
    end
    
    for trial = 1:length(trial_data)
        data = zeros(size(trial_data(trial).(signals{1,1}),1),sum(n_signals));
        count = 0;
        for i = 1:size(signals,1)
            temp_data = cat(1,trial_data(trial).(signals{i,1}));
            data(:,count+(1:n_signals(i))) = temp_data(:,signals{i,2});
            count = count + n_signals(i);
        end
        
        % project into low-D space and pick the num_dims requested
        temp_proj = (data - repmat(mu,size(data,1),1)) * w;
        if strcmpi(algorithm,'fa') && fa_orthogonalize % orthogonalize
            [temp_proj,~] = orthogonalize(temp_proj',w);
            temp_proj = temp_proj'; % that code uses dimensions as rows
        end
        if ~isempty(num_dims)
            trial_data(trial).([[sig_name{:}] '_' lower(algorithm)]) = temp_proj(:,1:num_dims);
        else
            trial_data(trial).([[sig_name{:}] '_' lower(algorithm)]) = temp_proj;
        end
    end
end