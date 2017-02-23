%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trial_data,glm_info] = getGLM(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
glm_name      =  'default';
in_signals    =  {};%{'name',idx; 'name',idx};
out_signals   =  {};%{'name',idx};
do_lasso      =  false;
lasso_lambda  =  0;
lasso_alpha   =  0;
train_idx     =  1:length(trial_data);
assignParams(who,params); % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if ~iscell(in_signals) || ~iscell(out_signals), error('input/output info must be in cells'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build inputs and outputs for training
x = get_vars(trial_data(train_idx),in_signals);
y = get_vars(trial_data(train_idx),out_signals);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit GLMs
b = zeros(size(x,2)+1,size(y,2));
for iVar = 1:size(y,2) % loop along outputs to predict
    if do_lasso % not quite implemented yet
        [b,s] = lassoglm(x,y(:,iVar),'poisson','lambda',lasso_lambda,'alpha',lasso_alpha);
        b = [s.Intercept; b];
    else
        [b(:,iVar),~,s(iVar)] = glmfit(x,y(:,iVar),'poisson');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add predictions to trial_data
for trial = 1:length(trial_data)
    x  = get_vars(trial_data(trial),in_signals);
    yfit = zeros(size(x,1),size(b,2));
    for iVar = 1:size(b,2)
        yfit(:,iVar) = exp([ones(size(x,1),1), x]*b(:,iVar));
    end
    trial_data(trial).(['glm_' glm_name]) = yfit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package up outputs
s = rmfield(s,{'resid','residp','residd','resida','wts'});
glm_info = struct( ...
    'glm_name',     glm_name, ...
    'in_signals',   {in_signals}, ...
    'out_signals',  {out_signals}, ...
    'b',            b, ...
    's',            s, ...
    'do_lasso',     do_lasso, ...
    'lasso_lambda', lasso_lambda, ...
    'lasso_alpha',  lasso_alpha);

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