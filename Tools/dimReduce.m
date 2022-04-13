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
%   'pca'     :  vanilla PCA, based on Matlab's built-in pca function
%   'ppca     :  probabalistc PCA, based on Matlab's built-in ppca function
%   'fa'      :  factor analysis, based on Matlab's factoran
%   'isomap'  :  nonlinear dimensionality reduction. See util/lib/isomap
%       NOTE: for all but PCA, must provide num_dims as input
%
% NOTE: centers data by default! Thus to reconstruct scores you need the
%       means of each signal. It stores these
%
% A hint: if you want the default PCA parameters, you can pass in just the
% signals that you want instead of a params struct.
%   e.g., trial_data = dimReduce(trial_data,'M1_spikes');
%
% INPUTS:
%   trial_data : the struct
%   params     : struct containing parameters
%     .algorithm      : (string) which algorithm, e.g., 'pca','fa','ppca'
%     .signals        : (cell) which signals. Two options:
%                           1) {'NAME1','NAME2',etc}
%                           2) {'NAME1',idx; 'NAME2',idx; etc}
%                                   Here idx is which columns to use
%                                   Note: can use 'all' as idx for all
%     .use_trials     : which trials to use for computing the low-D space (default: all)
%                         Note: when adding projections, does so for all
%                         trials in trial_data, not just use_trials
%     .num_dims       : how many dimensions (e.g. for PPCA, FA). For PCA,
%                        do nothing and it returns same dimensionality as
%                        input, or specify a value if you like.
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
% TO DO:
%   -  extend isomap to allow out of sample embedding
%
% Written by Matt Perich. Updated March 2019.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,info_out] = dimReduce(trial_data, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
algorithm       = 'pca';
signals         =  getTDfields(trial_data,'spikes');
use_trials      =  1:length(trial_data);
num_dims        =  []; % how many dimensions, needed for PPCA, FA, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters you can change that aren't described in header
sig_name          = '';     % output will be in field "SIG_NAMES_ALGORITHM". Defaults to concatenated names of signals
center_data       = true;   % whether to center data
sqrt_transform    = false;  % square root transform before reduction (projections don't have it)
do_smoothing      = false;  % will smooth before dim reduction  (trial_data projections are unsmoothed)
width             = 0.05;   %   gaussian kernel s.d. for smoothing
% PCA  parameters ---------------------------------------------------------
pca_algorithm     = 'svd';  % algorithm for PCA
pca_economy       = false;  % if num samples < degrees of freedom, will pad with zeros to keep output same size as degrees of freedom
% FA parameters -----------------------------------------------------------
fa_orthogonalize  = true;   % whether to orthogonalize the projections
fa_rotate         = 'none'; % rotation  to apply
private_var       = [];     % private variance for each channel (when fit by FA)
% Isomap parameters  ------------------------------------------------------
iso_n             = 12;     % neighborhood size (value for epsilon or k)
iso_function      = 'k';    % neighborhood function ('epsilon' or 'k')
iso_component     = 1;      % which connected component to embed, if more than one.
iso_display       = false;  % plot residual variance and 2-D embedding?
% parameters for projecting data ------------------------------------------
add_proj_to_td    = true;   % whether to add projections
recenter_for_proj = false;  % whether to recenter data before projecting into PC space
% initialize other parameters ---------------------------------------------
w                 = [];     % w is used to know if params was info_out (e.g. whether to recompute space)
mu                = [];     % mu is the mean from fitting, only filled if info_out is passed in
verbose           = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if ~isstruct(params)
        % assume the signals were provided
        signals = params;
    else
        if isfield(params,'trial_idx')
            warning('''trial_idx'' input has been changed to ''use_trials''... FYI.');
            use_trials = trial_idx;
        end
        assignParams(who,params); % overwrite parameters
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process and prepare inputs
trial_data  =  check_td_quality(trial_data);
signals     =  check_signals(trial_data(1),signals);
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
        td = smoothSignals(td,struct('signals',{signals},'width',width));
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
            [w, scores, eigen,~,~,mu] = pca(data,'Algorithm',pca_algorithm,'Centered',center_data,'Economy',pca_economy);
            stats = [];
        case 'ppca'
            if isempty(num_dims), error('Must specify number of dimensions when using PPCA'); end
            [w, scores, eigen, mu] = ppca(data,num_dims);
            stats = [];
        case 'fa'
            if isempty(num_dims), error('Must specify number of dimensions when using FA'); end
            [fa_params,~] = fastfa(data',num_dims);
            mu = fa_params.d';
            eigen = [];
            stats = [];
            private_var = fa_params.Ph';

            w = fa_params.L'; % note: this is the FA loading matrix, not a projection matrix from data to factor space
            [fa_post,~] = fastfa_estep(data',fa_params);
            scores = fa_post.mean';
            if fa_orthogonalize % orthogonalize
                [scores,w] = orthogonalize(scores',w');
                scores = scores'; % that code uses dimensions as rows
                w = w';
            end
        case 'isomap'
            if isempty(num_dims), error('Must specify number of dimensions when using Isomap'); end
            if length(use_trials) ~= length(trial_data)
                error('Sorry, for now isomap cannot project new data. The embedding must be calculated using all data.')
            end
            iso_opts = struct( ...
                'dims',1:num_dims, ...
                'comp',iso_component, ...
                'display',iso_display, ...
                'overlay',true,  ...
                'verbose',verbose);
            
            % center data
            if center_data
                mu = mean(data,1);
                data = data - repmat(mu,size(data,1),1);
            else
                mu = zeros(1,size(data,2));
            end
            
            % compute distances
            D = L2_distance(data', data');
            % find embedding
            
            if size(data,1) > 1000
                disp('------------------------------------')
                disp('Isomap may take a long time with > 1000 data points. Enabling verbosity so you can see.');
                iso_opts.verbose = true;
            end
            [Y, R, E, iso_info] = isomap(D, iso_function, iso_n, iso_opts);
            
            % package things up
            w = iso_info;
            scores = Y.coords{end}';
            stats.R = R;
            eigen = [];
            
            if length(Y.index) ~= size(data,1)
                disp('WARNING: Isomap could not find a single fully-connected component.');
            end
            
        otherwise
            error('Algorithm for dimReduce not recognized');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Package up outputs
    out_params = struct();
    out_params.signals     =  signals;
    out_params.use_trials  =  use_trials;
    info_out = struct('algorithm',algorithm,'w',w,'scores',scores,'eigen',eigen,'mu',mu,'private_var',private_var,'signals',{signals},'out_params',out_params,'sig_name',sig_name,'stats',stats,'num_dims',num_dims);
else
    info_out = params;
    if strcmpi(algorithm,'isomap')
        error('Isomap currently does not support embedding out of sample points.');
    end
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
    
    switch lower(algorithm)
        
        case 'isomap' % isomap has to go off the embedding for now
            scores = info_out.scores;
            count = 0;
            for trial = 1:length(trial_data)
                N = size(trial_data(trial).(signals{1,1}),1);
                
                % get the
                temp_proj = scores(count+1:count+N,:);
                
                % update counter
                count = count + N;
                
                if ~isempty(num_dims)
                    trial_data(trial).([[sig_name{:}] '_' lower(algorithm)]) = temp_proj(:,1:num_dims);
                else
                    trial_data(trial).([[sig_name{:}] '_' lower(algorithm)]) = temp_proj;
                end
            end

        case 'fa' % have to estimate the posterior from parameters
            for trial = 1:length(trial_data)
                % initialize data
                data = zeros(size(trial_data(trial).(signals{1,1}),1),sum(n_signals));
                
                % loop along all of the columns
                count = 0;
                for iSig = 1:size(signals,1)
                    temp_data = cat(1,trial_data(trial).(signals{iSig,1}));
                    data(:,count+(1:n_signals(iSig))) = temp_data(:,signals{iSig,2});
                    count = count + n_signals(iSig);
                end
                
                % estimate posterior from fa params
                [Z,~] = fastfa_estep(data',struct('L',w','Ph',private_var','d',mu'));
                trial_data(trial).([[sig_name{:}] '_' lower(algorithm)]) = Z.mean';
            end
            
        otherwise % project into the low-D space using the weight matrix
            if center_data
                if recenter_for_proj % find the new mean
                    mu = mean(get_vars(trial_data,signals),1);
                else % use the original mean
                    mu = info_out.mu;
                    if isempty(mu)
                        error('Could not find the original mean values for centering.');
                    end
                end
            else % don't subtract the mean
                mu = zeros(1,sum(n_signals));
            end
            
            for trial = 1:length(trial_data)
                % initialize data
                data = zeros(size(trial_data(trial).(signals{1,1}),1),sum(n_signals));
                
                % loop along all of the columns
                count = 0;
                for iSig = 1:size(signals,1)
                    temp_data = cat(1,trial_data(trial).(signals{iSig,1}));
                    data(:,count+(1:n_signals(iSig))) = temp_data(:,signals{iSig,2});
                    count = count + n_signals(iSig);
                end
                
                % project into low-D space and pick the num_dims requested
                temp_proj = (data - repmat(mu,size(data,1),1)) * w;
                if ~isempty(num_dims)
                    trial_data(trial).([[sig_name{:}] '_' lower(algorithm)]) = temp_proj(:,1:num_dims);
                else
                    trial_data(trial).([[sig_name{:}] '_' lower(algorithm)]) = temp_proj;
                end
            end
    end
end
