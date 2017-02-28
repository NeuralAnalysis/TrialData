%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data,pca_info] = getPCA(trial_data, params)
%
% [trial_data,pca_info] = getPCA(trial_data, params);
%   Computes PCA projection for time-varying data. Will add scores to each trial
% for use later. Must pass in 'signals' field for struct. Note that this can
% be a cell of strings to pool data from multiple arrays. In this case, the
% rows of w will be as if you concatenated the columns of signals together in the
% order they were provided.
%
% trial_data = getPCA(trial_data, pca_info);
%   Uses an old w from previous getPCA call to add scores to
% trial_data as SIGNALS{:}_pca.
%
% NOTE: centers data by default! Thus to reconstruct scores you need the
%       means of each signal.
%
% INPUTS:
%   trial_data : the struct
%   params     : struct containing parameters
%     .signals        : (cell) which signals. Two options:
%                           1) {'NAME1','NAME2',etc}
%                           2) {'NAME1',idx; 'NAME2',idx; etc}
%                                   Here idx is which columns to use
%                                   Note: can use 'all' as idx for all
%     .trial_idx      : which trials to use for PCA (default: all)
%                         Note: when adding projections, does so for all
%                         trials in trial_data, not just trial_idx
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
%   e.g. to add scores to trial_data later using the above output
%       trial_data = getPCA(trial_data, pca_info);
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,pca_info] = getPCA(trial_data, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
trial_idx       =  1:length(trial_data);
signals         =  getTDfields(trial_data,'spikes');
do_plot         =  false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters you can change that aren't described in header
do_smoothing    = false;  % will smooth before PCA  (trial_data projections are unsmoothed)
kernel_SD       = 0.05;  %   gaussian kernel s.d. for smoothing
pca_algorithm   = 'svd'; % algorithm for PCA
pca_centered    = true;  % whether to center data
add_proj_to_td  = true;  % whether to add PCA projections
w               = [];    % w is used to know if params was pca_info
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process and prepare inputs
signals = check_signals(trial_data(1),signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build PCA model for M1
if isempty(w)
    td = trial_data;
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Package up outputs
    pca_params = struct( ...
        'signals',{signals}, ...
        'trial_idx',trial_idx);
    pca_info = struct('w',w,'scores',scores,'eigen',eigen,'params',pca_params);
else
    pca_info = params;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add scores to trial_data
if add_proj_to_td
    sig_names = cellfun(@(x) strrep(x,'_spikes',''),signals(:,1),'uni',0);
    n_signals = cellfun(@(x) length(x),signals(:,2));
    
    if pca_centered
        mu = mean(get_vars(trial_data,signals),1);
    else
        mu = zeros(1,sum(n_signals));
    end
    
    for trial = 1:length(trial_data)
        data = zeros(size(trial_data(trial).(signals{1,1}),1),sum(n_signals));
        count = 0;
        for i = 1:size(signals,1)
            temp_data = cat(1,trial_data(trial).(signals{i,1}));
            data(:,count+(1:n_signals(i))) = temp_data(:,signals{i,2});
            count = count + n_signals(i);
        end
        
        trial_data(trial).([[sig_names{:}] '_pca']) = (data - repmat(mu,size(data,1),1)) * w;
    end
end