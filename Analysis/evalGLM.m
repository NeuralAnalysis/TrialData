%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%   trial_data : the struct
%   glm_name   : string name of glm to evaluate
%                   OR {'GLM1','GLM2'} to do relative pR2 of 1 rel to 2
function varargout = evalGLM(trial_data,glm_name,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
num_bootstraps   =  1000;
if nargin > 2, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if iscell(glm_name) % we are doing relative pr2
    if size(trial_data(1).([glm_name{1} '_glm']),2) ~= size(trial_data(1).([glm_name{2} '_glm']),2)
        error('Different numbers of variables for rpr2');
    end
    if length(glm_name) ~= 2
        error('For relative pR2, give glm_name as {''GLM1'',''GLM2''}');
    end
end

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
