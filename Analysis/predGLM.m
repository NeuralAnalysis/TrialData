%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = predGLM(trial_data,glm_info,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
num_bootstraps   =  1000;
rpr2_glm         =  []; % existing glm_name for rel pseudo-R2 calculation
do_single_trial  =  false;
if nargin > 2, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
glm_name    = glm_info.glm_name;
in_signals  = glm_info.in_signals;
out_signals = glm_info.out_signals;
b = glm_info.b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make predictions for each trial
if do_single_trial%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure out how many signals there will be
    for i = 1:size(out_signals,1)
        % if second entry is 'all', use all
        if ischar(out_signals{i,2})
            idx{i} = 1:size(td(1).(out_signals{i,1}),2);
        else
            idx{i} = out_signals{i,2};
        end
    end
    
    % Get performance by trial and add predictions to trial_data
    pr2 = NaN(length(trial_data),sum(cellfun(@(x) length(x),idx)),2);
    if ~isempty(rpr2_glm)
        rpr2 = NaN(length(trial_data),sum(cellfun(@(x) length(x),idx)),2);
    end
    
    for trial = 1:length(trial_data)
        x  = get_vars(trial_data(trial),in_signals);
        y = get_vars(trial_data(trial),out_signals);
        if ~isempty(rpr2_glm)
            yfit_rpr2 = get_vars(trial_data,{['glm_' rpr2_glm],'all'});
            [yfit,pr2(trial,:,:),rpr2(trial,:,:)] = eval_glm(b,x,y,yfit_rpr2,num_bootstraps);
        else
            [yfit,pr2(trial,:,:),~] = eval_glm(b,x,y,[],num_bootstraps);
        end
        
        trial_data(trial).(['glm_' glm_name]) = yfit;
    end
    
else%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get performance with all trials pooled together
    x  = get_vars(trial_data,in_signals);
    y  = get_vars(trial_data,out_signals);
    % if you want a relative pseudo R2, get the old predictions
    if ~isempty(rpr2_glm)
        yfit_rpr2 = get_vars(trial_data,{['glm_' rpr2_glm],'all'});
        [~,pr2,rpr2] = eval_glm(b,x,y,yfit_rpr2,num_bootstraps);
    else
        [~,pr2,~] = eval_glm(b,x,y,[],num_bootstraps);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package up outputs
varargout{1} = trial_data;
varargout{2} = pr2;
if ~isempty(rpr2_glm)
    varargout{3} = rpr2;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = get_vars(td,signals)
idx = cell(1,size(signals,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure out how many signals there will be
for i = 1:size(signals,1)
    % if second entry is 'all', use all
    if ischar(signals{i,2})
        idx{i} = 1:size(td(1).(signals{i,1}),2);
    else
        idx{i} = signals{i,2};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get datapoints
x = zeros(size(cat(1,td.pos),1),sum(cellfun(@(x) length(x),idx)));
count = 0;
for i = 1:size(signals,1)
    
    temp = cat(1,td.(signals{i,1}));
    x(:,count+(1:length(idx{i}))) = temp(:,idx{i});
    count = count + length(idx{i});
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yfit,pr2,rpr2] = eval_glm(b,x,y,yfit_rpr2,num_bootstraps)
yfit = zeros(size(y));
[pr2,rpr2] = deal(zeros(size(y,2),2));

% loop along output variables and compute
for iVar = 1:size(y,2)
    % predict from GLM output
    yfit(:,iVar) = exp([ones(size(x,1),1), x]*b(:,iVar));
    pr2(iVar,:) = get_pr2(y(:,iVar),yfit(:,iVar),num_bootstraps);
    
    if ~isempty(yfit_rpr2)
        rpr2(iVar,:) = get_rpr2(y(:,iVar),yfit_rpr2(:,iVar),yfit(:,iVar),num_bootstraps);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pr2 = get_pr2(y_test,y_fit,num_bootstraps)
% this is a really efficient way to bootstrap but you need temps
if num_bootstraps > 1
    bs = randi(length(y_test),length(y_test),num_bootstraps);
else
    bs = 1:size(y_fit,1);
end
pr2 = prctile(compute_pseudo_R2(y_test(bs),y_fit(bs),mean(y_test)),[2.5 97.5]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rpr2 = get_rpr2(y_test,y_fit_basic,y_fit_full,num_bootstraps)
% this is a really efficient way to bootstrap but you need temps
if num_bootstraps > 1
    bs = randi(length(y_test),length(y_test),num_bootstraps);
else
    bs = 1:size(y_fit,1);
end
rpr2 = prctile(compute_rel_pseudo_R2(y_test(bs),y_fit_basic(bs),y_fit_full(bs)),[2.5 97.5]);
end
