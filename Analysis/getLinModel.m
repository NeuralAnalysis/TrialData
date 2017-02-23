%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [orig,pred,c,b] = linLinModel(trial_data,params)
%
% 	Compute linear predictions of a signal. For example. predict X/Y
% position using M1 firing rates. Doesn't incorporate history outright. If
% you want history you should use dupeAndShift or convBasisFunc.
%
%   [trial_data,model_info] = getLinModel(trial_data,params);
%       This mode will fit a model and return the struct with predictions
%       added as a field, as well as parameters.
%
%   trial_data            = getLinModel(trial_data,model_info)
%       This mode will take the model_info from a previous getLinModel call and
%       add predictions as a field to trial_data (e.g. for new trials)
%
% WORK IN PROGRESS. To do:
%   1) Implement something like ridge regression for regularization
%   2) Merge with getGLM somehow, since they're very similar?
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%     .model_name  : (string) name for this model
%     .in_signals  : {'NAME1',idx1; 'NAME2', idx2}
%     .out_signals : {'NAME',idx}
%     .train_idx   : which trials for training (wll add pred to all)
%     .fit_metric  : which metric for evaluating fit ('corr','r2','vaf')
%
% OUTPUTS:
%   orig : time-varying original signal
%   pred : tim-varying predicted signal
%   c    : output of fit_metric
%   b    : weights computed for model
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,model_info] = getLinModel(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
model_name   =  '';
in_signals   =  {};
out_signals  =  {};
b            =  [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented parameters to overwrite
add_pred_to_td       =  true;
td_fieldname_prefix  =  'linmodel';
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if isempty(in_signals) || isempty(out_signals) || ~iscell(in_signals) || ~iscell(out_signals)
    error('input/output info must be provided in cells');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit model
if isempty(b)
    in_data = get_vars(trial_data,in_signals);
    out_data = get_vars(trial_data,out_signals);
    
    b = zeros(size(in_data,2),size(out_data,2));
    for iDim = 1:size(out_data,2)
        b(:,iDim) = [ones(size(in_data,1),1), in_data]\out_data(:,iDim);
    end
else % use an old model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fn = fieldnames(params);
    for i = 1:length(fn)
        assignin('caller',fn{i},params.(fn{i}));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add to trial_data struct
if add_pred_to_td
    for trial = 1:length(trial_data)
        in_data = get_vars(trial_data(trial),in_signals);
        yfit = zeros(size(in_data,1),size(b,2));
        for iVar = 1:size(b,2)
            yfit(:,iVar) = [ones(size(in_data,1),1), in_data]*b(:,iVar);
        end
        trial_data(trial).([td_fieldname_prefix '_' model_name]) = yfit;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package up outputs
model_info = struct( ...
    'model_name',   model_name, ...
    'in_signals',   {in_signals}, ...
    'out_signals',  {out_signals}, ...
    'train_idx',    train_idx, ...
    'b',            b);

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