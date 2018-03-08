%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trial_data,model_info] = getModel(trial_data,params)
%
%   This function will fit a linear model to take any combination of inputs
% and predict some outputs. Supports GLM or normal linear models. Doesn't
% incorporate history outright. If you want history you should use
% dupeAndShift or convBasisFunc.
%
% Actually has two modes:
%
%   [trial_data,model_info] = getModel(trial_data,model_type,params);
%       This mode will fit a model and return the struct with GLM predictions
%       added as a field, as well as GLM parameters.
%
%   trial_data            = getModel(trial_data,model_type,model_info)
%       This mode will take the model_info from a previous getModel call and
%       add model predictions as a field to trial_data (e.g. for new trials)
%
%
% INPUTS:
%   trial_data : the struct
%   in_struct  : a struct of inputs. Can be one of two things:
%       1) glm_info   : struct of model fit info from getGLM call (for predicting)
%       2) params     : parameter struct (NOTE: must NOT have 'b' or 's' fields)
%            .model_type   : (string) which model to use (REQUIRED)
%                               'linmodel' : a simple linear filter
%                               'glm'      : a generalized linear model
%            .model_name   : (string) unique name for this model fit
%            .in_signals   : (cell) GLM inputs in form {'name',idx; 'name',idx};
%            .out_signals  : (cell) GLM outputs in form {'name',idx}
%            .train_idx    : trial indices for training (will predict for all)
%          GLM-specific parameters:
%            .do_lasso     : flag to use lasso regularization (default: false)
%            .lasso_lambda : lambda parameter for lasso
%            .lasso_alpha  : alpha parameter for lasso
%
% Note: at the moment, only supports one output field. But can have any
%       arbitrary number of input fields, which it will stick together
%
% To do:
%   1) Implement something like ridge regression for linmodel
%
% EXAMPLES:
%   e.g. to fit a GLM
%       [trial_data,glm_info] = getGLM(trial_data(train_trials),params);
%   e.g. to use a previous model make predictions for new data
%       trial_data            = getGLM(trial_data(new_trials),glm_info)
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,model_info] = getModel(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
model_type    =  '';
model_name    =  'default';
in_signals    =  {};%{'name',idx; 'name',idx};
out_signals   =  {};%{'name',idx};
train_idx     =  1:length(trial_data);
% GLM-specific parameters
do_lasso      =  false;
lasso_lambda  =  0;
lasso_alpha   =  0;
nn_params = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here are some parameters that you can overwrite that aren't documented
add_pred_to_td       =  true;        % whether to add predictions to trial_data
glm_distribution     =  'poisson';   % which distribution to assume for GLM
td_fn_prefix         =  '';  % name prefix for trial_data field
b                    =  [];          % b and s identify if model_info was
s                    =  [];          %    provided as a params input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if isempty(model_type), error('Must specify what type of model to fit'); end
if isempty(in_signals) || isempty(out_signals)
    error('input/output info must be provided');
end
if isempty(td_fn_prefix), td_fn_prefix = model_type; end
in_signals = check_signals(trial_data(1),in_signals);
out_signals = check_signals(trial_data(1),out_signals);
if iscell(train_idx) % likely to be meta info
    train_idx = getTDidx(trial_data,train_idx{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(b)  % fit a new model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % build inputs and outputs for training
    x = get_vars(trial_data(train_idx),in_signals);
    y = get_vars(trial_data(train_idx),out_signals);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit GLMs
    b = zeros(size(x,2)+1,size(y,2));
    for iVar = 1:size(y,2) % loop along outputs to predict
        switch lower(model_type)
            case 'glm'
                if do_lasso % not quite implemented yet
                    % NOTE: Z-scores here!
                    [b_temp,s_temp] = lassoglm(zscore(x),y(:,iVar),glm_distribution,'lambda',lasso_lambda,'alpha',lasso_alpha);
                    b(:,iVar) = [s_temp.Intercept; b_temp];
                else
                    [b(:,iVar),~,s_temp] = glmfit(x,y(:,iVar),glm_distribution);
                end
                if isempty(s)
                    s = s_temp;
                else
                    s(iVar) = s_temp;
                end
            case 'linmodel'
                b(:,iVar) = [ones(size(x,1),1), x]\y(:,iVar);
            case 'neural_net'
                trained_net = trainNetwork(x, y, nn_params,  
            end
    end
else % use an old GLM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % these parameters should already be assigned from assignParams
%     fn = fieldnames(params);
%     for i = 1:length(fn)
%         assignin('caller',fn{i},params.(fn{i}));
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add predictions to trial_data
if add_pred_to_td
    for trial = 1:length(trial_data)
        x  = get_vars(trial_data(trial),in_signals);
        
        yfit = zeros(size(x,1),size(b,2));
        for iVar = 1:size(b,2)
            switch lower(model_type)
                case 'glm'
                    if do_lasso
                        yfit(:,iVar) = exp([ones(size(x,1),1), zscore(x)]*b(:,iVar));
                    else
                        yfit(:,iVar) = exp([ones(size(x,1),1), x]*b(:,iVar));
                    end
                case 'linmodel'
                    yfit(:,iVar) = [ones(size(x,1),1), x]*b(:,iVar);
            end
        end
        trial_data(trial).([td_fn_prefix '_' model_name]) = yfit;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package up outputs
switch lower(model_type)
    case 'glm'
        if strcmpi(model_type,'glm') && ~do_lasso
            if isfield(s,'resid')
                s = rmfield(s,{'resid','residp','residd','resida','wts'});
            end
        end
        model_info = struct( ...
            'model_type',   model_type, ...
            'model_name',   model_name, ...
            'in_signals',   {in_signals}, ...
            'out_signals',  {out_signals}, ...
            'train_idx',    train_idx, ...
            'b',            b, ...
            's',            s, ...
            'do_lasso',     do_lasso, ...
            'lasso_lambda', lasso_lambda, ...
            'lasso_alpha',  lasso_alpha);
        
    case 'linmodel'
        model_info = struct( ...
            'model_type',   model_type, ...
            'model_name',   model_name, ...
            'in_signals',   {in_signals}, ...
            'out_signals',  {out_signals}, ...
            'train_idx',    train_idx, ...
            'b',            b);
end