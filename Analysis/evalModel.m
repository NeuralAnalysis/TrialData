%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function metric = evalModel(trial_data,params)
%
%   Evaluates quality of model fit. Input trial_data must have
% recognizable field added to use for the predictions. You need to have the
% model_type, model_name, and out_signals parameters defined at input, as
% well as the eval_metric input.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .out_signals  : which output signal is predicted by model_name
%       .model_name   : string name of model to evaluate
%                         OR {'NAME1','NAME2'} to do relative metric of 2 rel to 1
%       .trial_idx    : trials to evaluate. Ways to use:
%                     1) 1:end treats each trial separately
%                     2) 1:N:end predicts in bins of size N trials
%                     3) [1,end] returns a single value for predicting all trials
%                         DEFAULT: [1,length(trial_data]
%       .eval_metric  : (string) name of metric for evaluation
%                           'pr2' : pseudo-R2 (Default)
%                           'vaf' : VAF (good thing I have these descriptions)
%                           'r2'  : R2 (as in, square of r)
%                           'r'   : correlation coefficient
%       .block_trials : if true, takes input of trial indices and pools
%                       them together for a single eval. If false, treats the trial indices
%                       like a list of blocked testing segments
%       .num_boots    : # bootstrap iterations to use (if <2, doesn't bootstrap)
%
% OUTPUTS:
%   metric : calculated evaluation metric
%                Note: will return relative metric if model_name is 1x2 cell of names
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metric = evalModel(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
model_type       =  '';
out_signals      =  [];
model_name       =  [];
trial_idx        =  [1,length(trial_data)];
eval_metric      =  '';
block_trials     =  false;
num_boots        =  1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented parameters
td_fn_prefix     =  '';    % prefix for fieldname
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
possible_metrics = {'pr2','vaf','r2','r'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if isempty(td_fn_prefix), td_fn_prefix = model_type; end
if isempty(model_type), error('Unknown model type.'); end
if isempty(out_signals), error('Need to provide output signal'); end
if isempty(model_name), error('No model_name provided'); end
if isempty(eval_metric), error('Must provide evaluation metric.'); end
if iscell(model_name) && length(model_name) == 2 % we are doing relative metric
    if size(trial_data(1).([td_fn_prefix '_' model_name{1}]),2) ~= size(trial_data(1).([td_fn_prefix '_' model_name{2}]),2)
        error('Different numbers of variables for relative metric calc');
    end
    if length(model_name) ~= 2
        error('For relative metric, give model_name as {''NAME1'',''NAME2''}');
    end
end
if ~iscell(model_name), model_name = {model_name}; end
if ~any(ismember(eval_metric,possible_metrics)), error('Metric not recognized.'); end
out_signals = check_signals(trial_data(1),out_signals);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate pr2
if block_trials
    metric = NaN(sum(cellfun(@(x) length(x),out_signals(:,2))),2);
    y = get_vars(trial_data(trial_idx),out_signals);
    yhat1 = cat(1,trial_data(trial_idx).([td_fn_prefix '_' model_name{1}]));
    if length(model_name) == 1
        for iVar = 1:size(y,2)
            metric(iVar,:) = get_metric(y(:,iVar),yhat1(:,iVar),eval_metric,num_boots);
        end
    else % relative metric
        yhat2 = cat(1,trial_data(trial_idx).([td_fn_prefix '_' model_name{2}]));
        for iVar = 1:size(y,2)
            metric(iVar,:) = get_metric(y(:,iVar),yhat1(:,iVar),yhat2(:,iVar),eval_metric,num_boots);
        end
    end
    if num_boots < 2
        metric = metric(:,1);
    end
else
    % quick hack here for when we aren't binning and want single trials
    if unique(diff(trial_idx)) == 1, trial_idx = [trial_idx, trial_idx(end)+1]; end
    metric = NaN(length(trial_idx)-1,sum(cellfun(@(x) length(x),out_signals(:,2))),2);
    for i = 1:length(trial_idx)-1
        trials = trial_idx(i):trial_idx(i+1)-1;
        
        y = get_vars(trial_data(trials),out_signals);
        yhat1 = cat(1,trial_data(trials).([td_fn_prefix '_' model_name{1}]));
        if length(model_name) == 1
            for iVar = 1:size(y,2)
                metric(i,iVar,:) = get_metric(y(:,iVar),yhat1(:,iVar),eval_metric,num_boots);
            end
        else % relative metric
            yhat2 = cat(1,trial_data(trials).([td_fn_prefix '_' model_name{2}]));
            for iVar = 1:size(y,2)
                metric(i,iVar,:) = get_metric(y(:,iVar),yhat1(:,iVar),yhat2(:,iVar),eval_metric,num_boots);
            end
        end
    end
    % if we don't need 3-D, ditch the single dim
    if num_boots < 2
        metric = squeeze(metric(:,:,1));
    end
end



end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metric = get_metric(varargin)
y_test         = varargin{1};
eval_metric    = varargin{end-1};
num_bootstraps = varargin{end};

% this is a really efficient way to bootstrap...
if num_bootstraps > 1
    bs = randi(length(y_test),length(y_test),num_bootstraps);
else
    bs = 1:length(y_test);
end

y_fit = varargin{2};
if nargin == 4 % normal metric
    switch lower(eval_metric)
        case 'pr2'
            metric = compute_pseudo_R2(y_test(bs),y_fit(bs),mean(y_test));
        case 'vaf'
            metric = compute_vaf(y_test(bs),y_fit(bs));
        case 'r2'
            metric = compute_r2(y_test(bs),y_fit(bs));
    end
elseif nargin == 5 % relative metric
    y_fit2 = varargin{3};
    switch lower(eval_metric)
        case 'pr2'
            metric = compute_rel_pseudo_R2(y_test(bs),y_fit(bs),y_fit2(bs));
        case 'vaf'
            error('VAF not yet implemented for relative metrics');
        case 'r2'
            error('R2 not yet implemented for relative metrics');
        case 'r'
            error('r not yet implemented for relative metrics');
    end
end
% no need to repeat
if num_bootstraps > 1
    metric = prctile(metric,[2.5,97.5]);
end
end