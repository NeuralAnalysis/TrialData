%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  example script to show how to use TrialData  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TD structs can have any arbitrary fields, which contain whatever you
% want. However, you should respect these constraints:
%   1) All time-varying signals must be on a unified time vector.
%   2) try not to have nested structs. The code is built to be able to
%   rapidly/easily query and partition the data because it's all accessible
%   within the top level of the struct
%
% With this, analysis becomes extremely easy. The framework is very
% flexible, and all of the "standard" functions (binning, smoothing,
% trimming, aligning) are built-in.

%  load in the struct
load('ex_td.mat');
% this struct contains Chewie M1/PMd recordings
%   It's a single recording with all of the  signals (kinematics, etc)
%   synchronized, and event indices corresponding to the bins for target on
%   etc. There is also associated metadata (you can add a field with
%   whatever arbitrary data you want).

% the standard processing tools have standardized inputs/outputs.
%   trial_data = FUNCTION(trial_data,params)
%
% This is one of the key features of the framework... every function (with
% few exceptions) takes the same input and spits out the same struct. So
% it's completely modular.
%
% params is a struct where each field is a parameter name. The names
% have to match whatever is in the function... see the function headers.
% But you can make one struct and pass it around to many functions, there
% isn't harm to having unnecessary fields. This is convenient for storing
% parameters with results/data.
%
% Some functions are pre-programmed to give you access to one parameters
% instead of a struct. You just pass in the value of that parameter; other
% parameters are set to defaults. Or you can pass in a struct and overwrite
% anything.

% the structs tend to be memory efficient so we can keep the original
% and do our destructive edits on a copy.
td = trial_data;

% you can first perform any processing you want on the full dataset
%   e.g. smoothing... things that might have edge effects

% we want to compute smooth firing rates, so convolve spikes with a
% gaussian kernel. If you like the default value of 50 ms, just tell it the
% names of the signals
% td = smoothSignals(td,{'M1_spikes','PMd_spikes'});
% maybe you want 100ms width instead of 50, so you can specify:
td = smoothSignals(td,struct('signals',{{'M1_spikes','PMd_spikes'}},'width',0.1));

% we started at 10ms bins, make it 20ms
td = binTD(td,2);

% now we want to break it up into trials. This is easily done based on the
% events. For my data, idx_trial_start is useful for this.
% Often there is metadata associated with trials, so  you have to manually
% tell it that these fields are linked. idx_ fields are naturally split.
% They get filled with NaN if none are found in the range that defines the
% new trial.
td = splitTD(td,struct( ...
    'split_idx_name','idx_trial_start', ...
    'linked_fields',{{'result','target_direction'}}));

% but we don't want every trial. Only get the rewarded ones
[~,td] = getTDidx(td,'result','R');
% this function breaks convention a bit, because the first output is
% indices meeting the queries. This allows inline indexing of various
% arrays/functions, and if you just want the trials from the struct you
% only ask for the second output, as above


% sometimes the data gets corrupted or something and the limblab recording
% setup misses a word or whatever. Then the idx_ could be padded with NaN.
% You can easily remove these trials.
td = removeBadTrials(td);
% that function has more features... but the basic (no inputs required) is
% just to clean out any trials with NaN in idx fields.


% okay, let's say we only want to consider the task-relevant window. By
% simply splitting, it starts at one trial and ends at the bin before
% the next one. But that includes the return movement, which we aren't
% studying. So, trim to a more reasonable window. Let's say from 5 bins
% before target onset until the reward
td = trimTD(td, {'idx_target_on',-5},'idx_trial_end');
% note you can use things like 'start' and 'end', or cite specific idx_


% for later analyses, we don't want to low-firing neurons. You can easily
% remove them as below. The function can do other things (shunting etc)
td = removeBadNeurons(td,struct('min_fr',0.5));

% now we can stash this pre-processed struct so we can do further
% destructive edits
td_store = td;



%% now it's processed, let's have fun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example #1: compare preparation and movement subspaces
arrays = {'M1','PMd'};
pca_params = struct( ...
    'signals',{cellfun(@(x) [x '_spikes'],arrays,'uni',0)});

% first get planning activity
% dimReduce does dimensionality reduction... currently supports the basics
%   it returns a TD struct with an added field containing the projections
%   from the PCA (for example), and as a second output parameters related
%   to that (e.g. weight matrices)
[~,pca_prep] = dimReduce(trimTD(td_store,'idx_target_on','idx_go_cue'),pca_params);
% also note that because the outputs of every function are the same struct
% that you put into another function, you can stack them up inline like
% the above line. If you're really into compact code. I think doing it too
% much gets confusing for others though...

% now get movement activity subspace using trimTD into dimReduce
[~,pca_move] = dimReduce(trimTD(td_store,'idx_go_cue','end'),pca_params);

% now we can use weights defined within each epoch, and project the data
% from the full trial onto the axes we found in the above analysis.
% dimReduce does double duty here. If you pass in a weight matrix it will
% just apply it to the signals for all trials instead of computing a new
% one.
td_td_prep   = dimReduce(td_store,setfield(pca_params,'w',pca_prep.w));
td_td_move   = dimReduce(td_store,setfield(pca_params,'w',pca_move.w));


% for plotting, we could align on go cue with the following
% td_td_prep = trimTD(td_td_prep,{'idx_go_cue',-25},{'idx_go_cue',25});
% td_td_move = trimTD(td_td_move,{'idx_go_cue',-25},{'idx_go_cue',25});


% now make a plot
figure;
subplot(2,1,1); hold all;
for i = 1:length(td_td_prep)
    temp = td_td_prep(i).([[arrays{:}] '_pca']);
    plot(temp(:,1),'LineWidth',1,'Color',0.5*[1 1 1]);
    plot(td_td_prep(i).idx_go_cue,temp(td_td_prep(i).idx_go_cue,1),'o','Color','r');
end
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14);
ylabel('First PC Proj','FontSize',14)
title('Preparatory Subspace','FontSize',14);
subplot(2,1,2); hold all;
for i = 1:length(td_td_move)
    temp = td_td_move(i).([[arrays{:}] '_pca']);
    plot(temp(:,1),'LineWidth',1,'Color',0.5*[1 1 1])
    plot(td_td_move(i).idx_go_cue,temp(td_td_move(i).idx_go_cue,1),'o','Color','r');
end
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14);
ylabel('First PC Proj','FontSize',14);
title('Movement Subspace','FontSize',14);

% you should see a jumbled mess of trials since they aren't aligned on
% movement, but the movement subspace should
% be inactive in the left of the plot (more preparation)


%% next example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example #2: decoding!
td = td_store;

% let's predict velocity and speed. Velocity is already there
td = getNorm(td,struct('signals','vel','norm_name','speed'));

% we want a Wiener cascade with history. dupeAndShift adds a _shift field
% that contains different leads/lags of a signal
td = dupeAndShift(td,'M1_spikes',-(1:10),'PMd_spikes',-(1:10));

% we want a testing and training set for cross validation. This can simply
% be done by indexing trials
test_idx = 1:10;
train_idx = 11:length(td);

% getModel will build the wiener cascade. prepare the inputs
model_params = struct( ...
    'in_signals',{{'M1_spikes','PMd_spikes'}}, ...
    'out_signals',{{'vel','speed'}}, ...
    'model_type','linmodel', ...
    'model_name','vel', ...
    'train_idx',train_idx, ...
    'polynomial',3);
%   Note that I give it the training idx. It will build the model using
%   only those trials, even though I pass in the struct with all trials.
%   But the cool thing is that it then can add a field to the TD using all
%   trials which contains the prediction output of the model. So, by
%   indexing at test_idx, you instantly have your cross-validated
%   predictions
td = getModel(td,model_params);


% quantify how well it did. There is an evalModel function which has some
% stuff built in, but it's also a work in progress. We can access any
% signal we want directly using getSig
x_vel = getSig(td(test_idx),{'vel',1});
y_vel = getSig(td(test_idx),{'vel',2});
y_sp = getSig(td(test_idx),'speed');
x_vel_pred = getSig(td(test_idx),{'linmodel_vel',1});
y_vel_pred = getSig(td(test_idx),{'linmodel_vel',2});
y_sp_pred = getSig(td(test_idx),{'linmodel_vel',3});

vaf_x = compute_vaf(x_vel,x_vel_pred);
vaf_y = compute_vaf(y_vel,y_vel_pred);
vaf_s = compute_vaf(y_sp,y_sp_pred);

% now plot
figure;
ax(1) = subplot(3,1,1); hold all;
plot(x_vel,'LineWidth',2);
plot(x_vel_pred,'LineWidth',2);
title(['VAF = ' num2str(vaf_x,3)]);

ax(2) = subplot(3,1,2); hold all;
plot(y_vel,'LineWidth',2);
plot(y_vel_pred,'LineWidth',2);
title(['VAF = ' num2str(vaf_y,3)]);

ax(3) = subplot(3,1,3); hold all;
plot(y_sp,'LineWidth',2);
plot(y_sp_pred,'LineWidth',2);
title(['VAF = ' num2str(vaf_s,3)]);

h = legend({'Actual','Predicted'},'Location','SouthEast');
set(h,'Box','off');

linkaxes(ax,'x');




