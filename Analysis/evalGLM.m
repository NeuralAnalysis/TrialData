%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pr2 = evalGLM(trial_data,params)
%
%   Evaluates quality of GLM fit with pseudo-R2. Input trial_data must have
% 'glm_GLM_NAME' field added to use for the predictions. Params input must
% contain an out_signals so that the code knows what to compare the fit to.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .out_signals : which output signal is predicted by glm_name
%       .glm_name    : string name of glm to evaluate
%                         OR {'GLM1','GLM2'} to do relative pR2 of 2 rel to 1
%                         Default is to find the _glm field (error if more than one)
%       .trial_idx   : trials to evaluate. Ways to use:
%                     1) 1:end treats each trial separately
%                     2) 1:N:end predicts in bins of size N trials
%                     3) [1,end] returns a single value for predicting all trials
%                         DEFAULT: [1,length(trial_data]
%       .eval_metric : (string) name of metric for evaluation
%                           'pr2' : pseudo-R2 (Default)
%                           'vaf' : VAF
%                           'r2'  : R2 (good thing I have these descriptions)
%       .num_boots   : # bootstrap iterations to use (if <2, doesn't bootstrap)
%
% OUTPUTS:
%   pr2 : calculated evaluation metric (only psuedo-R2 for now)
%           Note: will return relative metric if glm_name is 1x2 cell of names
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pr2 = evalGLM(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
out_signals      =  [];
glm_name         =  [];
trial_idx        =  [1,length(trial_data)];
eval_metric      =  'pr2';
num_boots        =  1000;
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if isempty(out_signals), error('Need to provide output signal'); end
if ~iscell(out_signals), out_signals = {out_signals}; end
if iscell(glm_name) % we are doing relative pr2
    if size(trial_data(1).(['glm_' glm_name{1}]),2) ~= size(trial_data(1).(['glm_' glm_name{2}]),2)
        error('Different numbers of variables for rpr2');
    end
    if length(glm_name) ~= 2
        error('For relative pR2, give glm_name as {''GLM1'',''GLM2''}');
    end
end
if isempty(glm_name)
    fn = fieldnames(trial_data);
    idx = ~cellfun(@isempty,strfind(fn,'glm_'));
    if sum(idx) > 1, error('Multiple GLM models found. Please specify which to use.'); end
    glm_name = fn(idx);
end
if ~iscell(glm_name), glm_name = {glm_name}; end
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
pr2 = NaN(length(trial_idx)-1,sum(cellfun(@(x) length(x),idx)),2);
for i = 1:length(trial_idx)-1
    trials = trial_idx(i):trial_idx(i+1)-1;
    
    temp = cat(1,trial_data(trials).(out_signals{1}));
    temp1 = cat(1,trial_data(trials).(['glm_' glm_name{1}]));
    if length(glm_name) == 1 % pr2
        for iVar = 1:size(temp,2)
            pr2(i,iVar,:) = get_pr2(temp(:,iVar),temp1(:,iVar),eval_metric,num_boots);
        end
    else % relative pr2
        temp2 = cat(1,trial_data(trials).(['glm_' glm_name{2}]));
        for iVar = 1:size(temp,2)
            pr2(i,iVar,:) = get_pr2(temp(:,iVar),temp1(:,iVar),temp2(:,iVar),eval_metric,num_boots);
        end
    end
end

% if we don't need 3-D, ditch the single dim
pr2 = squeeze(pr2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pr2 = get_pr2(varargin)
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
            pr2 = prctile(compute_pseudo_R2(y_test(bs),y_fit(bs),mean(y_test)),[2.5 97.5]);
        case 'vaf'
            pr2 = prctile(compute_vaf(y_test(bs),y_fit(bs),mean(y_test)),[2.5 97.5]);
        case 'r2'
            pr2 = prctile(compute_r2(y_test(bs),y_fit(bs),mean(y_test)),[2.5 97.5]);
        otherwise, error('eval_metric not recognized.');
    end
elseif nargin == 5 % relative metric
    y_fit2 = varargin{3};
    switch lower(eval_metric)
        case 'pr2'
            pr2 = prctile(compute_rel_pseudo_R2(y_test(bs),y_fit(bs),y_fit2(bs)),[2.5 97.5]);
        case 'vaf'
            error('VAF not yet implemented for relative metrics');
        case 'r2'
            error('R2 not yet implemented for relative metrics');
        otherwise, error('eval_metric not recognized.');
    end
end
end

