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
%       .out_signals : which output signal is predicted by model_name
%       .model_name  : string name of model to evaluate
%                         OR {'NAME1','NAME2'} to do relative metric of 2 rel to 1
%       .trial_idx   : trials to evaluate. Ways to use:
%                     1) 1:end treats each trial separately
%                     2) 1:N:end predicts in bins of size N trials
%                     3) [1,end] returns a single value for predicting all trials
%                         DEFAULT: [1,length(trial_data]
%       .eval_metric : (string) name of metric for evaluation
%                           'pr2' : pseudo-R2 (Default)
%                           'vaf' : VAF (good thing I have these descriptions)
%                           'r2'  : R2 (as in, square of r)
%                           'r'   : correlation coefficient
%       .num_boots   : # bootstrap iterations to use (if <2, doesn't bootstrap)
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
num_boots        =  1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented parameters
td_fn_prefix     =  ''; % prefix for fieldname
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
if iscell(model_name) % we are doing relative metric
    if size(trial_data(1).([td_fn_prefix '_' model_name{1}]),2) ~= size(trial_data(1).([td_fn_prefix '_' model_name{2}]),2)
        error('Different numbers of variables for relative metric calc');
    end
    if length(model_name) ~= 2
        error('For relative metric, give model_name as {''NAME1'',''NAME2''}');
    end
end
if ~iscell(model_name), model_name = {model_name}; end
if ~iscell(out_signals), out_signals = {out_signals}; end
if ~any(ismember(eval_metric,possible_metrics)), error('Metric not recognized.'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure out how many signals there will be
for i = 1:size(out_signals,1)
    % if second entry is 'all', use all
    if ischar(out_signals{i,2})
        idx{i} = 1:size(trial_data(1).(out_signals{i,1}),2);
    else
        idx{i} = out_signals{i,2};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate pr2
% quick hack here for when we aren't binning and want single trials
if unique(diff(trial_idx)) == 1, trial_idx = [trial_idx, trial_idx(end)+1]; end
metric = NaN(length(trial_idx)-1,sum(cellfun(@(x) length(x),idx)),2);
for i = 1:length(trial_idx)-1
    trials = trial_idx(i):trial_idx(i+1)-1;
    
    temp = get_vars(trial_data(trials),out_signals);
    temp1 = cat(1,trial_data(trials).([td_fn_prefix '_' model_name{1}]));
    if length(model_name) == 1
        for iVar = 1:size(temp,2)
            metric(i,iVar,:) = get_metric(temp(:,iVar),temp1(:,iVar),eval_metric,num_boots);
        end
    else % relative metric
        temp2 = cat(1,trial_data(trials).([td_fn_prefix '_' model_name{2}]));
        for iVar = 1:size(temp,2)
            metric(i,iVar,:) = get_metric(temp(:,iVar),temp1(:,iVar),temp2(:,iVar),eval_metric,num_boots);
        end
    end
end

% if we don't need 3-D, ditch the single dim
metric = squeeze(metric);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metric = get_metric(varargin)
y_test         = varargin{1};
eval_metric    = varargin{end-1};
num_bootstraps = varargin{end};

% this is a really efficient way to bootstrap but you need temps
if num_bootstraps > 1
    bs = randi(length(y_test),length(y_test),num_bootstraps);
else
    bs = 1:length(y_test);
end

y_fit = varargin{2};
if nargin == 4 % normal metric
    switch lower(eval_metric)
        case 'pr2'
            metric = prctile(compute_pseudo_R2(y_test(bs),y_fit(bs),mean(y_test)),[2.5 97.5]);
        case 'vaf'
            metric = prctile(compute_vaf(y_test(bs),y_fit(bs)),[2.5 97.5]);
        case 'r2'
            metric = prctile(compute_r2(y_test(bs),y_fit(bs)),[2.5 97.5]);
        case 'r'
            metric = sqrt(prctile(compute_r2(y_test(bs),y_fit(bs)),[2.5 97.5]));
    end
elseif nargin == 5 % relative metric
    y_fit2 = varargin{3};
    switch lower(eval_metric)
        case 'pr2'
            metric = prctile(compute_rel_pseudo_R2(y_test(bs),y_fit(bs),y_fit2(bs)),[2.5 97.5]);
        case 'vaf'
            error('VAF not yet implemented for relative metrics');
        case 'r2'
            error('R2 not yet implemented for relative metrics');
        case 'r'
            error('R2 not yet implemented for relative metrics');
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = get_vars(td,signals)
idx = cell(1,size(signals,1));

% figure out how many signals there will be
for i = 1:size(signals,1)
    % if second entry is 'all', use all
    if ischar(signals{i,2})
        idx{i} = 1:size(td(1).(signals{i,1}),2);
    else
        idx{i} = signals{i,2};
    end
end

% get datapoints
x = zeros(size(cat(1,td.pos),1),sum(cellfun(@(x) length(x),idx)));
count = 0;
for i = 1:size(signals,1)
    
    temp = cat(1,td.(signals{i,1}));
    x(:,count+(1:length(idx{i}))) = temp(:,idx{i});
    count = count + length(idx{i});
end

end
