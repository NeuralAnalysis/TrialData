% %% Test building GPFA model
% root_dir = 'F:\trial_data_files\';
% 
% params.bin_w = 30;
% params.xdim = 8;
% params.kernsd = 30;
% params.arrays = 'M1';
% 
% dataSummary;
% 
% for iFile = 1:size(sessionList,1)
%     filename = [sessionList{iFile,1} '_CO_VR_' sessionList{iFile,2}];
%     if exist(fullfile(root_dir,[filename '_gpfa.mat']),'file')
%         load(fullfile(root_dir,[filename '_gpfa.mat']));
%     else
%         load(fullfile(root_dir,[filename '.mat']));
%     end
%     
%     % run the GPFA processing
%     [trial_data, gpfa_results] = run_gpfa(trial_data,params);
%     
%     % store the models for later use
%     gpfa_models.([params.arrays '_gpfa']) = gpfa_results.model;
%     gpfa_params.([params.arrays '_gpfa']) = gpfa_results.params;
%     
%     % save in a new file
%     save(fullfile(root_dir,[filename '_gpfa.mat']),'trial_data','gpfa_models','gpfa_params');
% end
% 
% %% Load in the just created gpfa files and compile to a master file
% all_trial_data = [];
% for iFile = 1:size(sessionList,1)
%     load(fullfile(root_dir,[sessionList{iFile,1} '_CO_VR_' sessionList{iFile,2} '_gpfa.mat']),'trial_data');
%     
%     if ~isfield(trial_data,'M1_spikes');
%         [trial_data.M1_spikes] = deal([]);
%     end
%     if ~isfield(trial_data,'PMd_spikes');
%         [trial_data.PMd_spikes] = deal([]);
%     end
%     if ~isfield(trial_data,'M1_gpfa');
%         [trial_data.M1_gpfa] = deal([]);
%     end
%     if ~isfield(trial_data,'PMd_gpfa');
%         [trial_data.PMd_gpfa] = deal([]);
%     end
%     
%     all_trial_data = [all_trial_data trial_data];
% end
% 
% clear trial_data;
% trial_data = all_trial_data;
% save(fullfile(root_dir,'All_CO_VR_gpfa.mat'),'trial_data')

%% Now do some analysis and plotting
clc;
close all;

epochs = {'BL','AD','WO'};

% gpfa parameters
data_bin_w = 10;
params.bin_w = 30;
params.xdim = 8;
params.kernsd = 30;
use_array = 'M1';

% identify a window around some event
% align_event = 'idx_go_cue';
% align_window = [0 0.3];
align_event = 'idx_movement_on';
align_window = [-0.4 0.2];
% align_event = 'idx_movement_on';
% align_window = [0 0.21];
% align_event = 'idx_reward';
% align_window = [-1.1 -0.5];
% align_event = 'idx_target_on';
% align_window = [0 0.6];

num_err_bins = 10;

dataSummary;
use_dates = sessionList(:,2);

results = repmat(struct(),1,length(use_dates));
for iDate = 1:length(use_dates)
    % get all trials for a given session
    date_trial_inds = find(get_trial_data_indices(trial_data,'date',use_dates{iDate}));
    td = trial_data(date_trial_inds);
    
    % get list of target directions
    utheta = unique([td.target_direction]);
    
    %%%%%%%%%%%%%%%%
    % Do Behavioral adaptation stuff
    
    % get average behavioral error for each trial in baseline
    dir_errs = zeros(1,length(utheta));
    for iDir = 1:length(utheta)
        idx = find(get_trial_data_indices(td,'target_direction',utheta(iDir),'epoch','BL'));
        
        % loop along trials
        all_errs = zeros(1,length(idx));
        for iTrial = 1:length(idx)
            temp = td(idx(iTrial)).pos;
            temp = temp(td(idx(iTrial)).idx_movement_on:td(idx(iTrial)).idx_movement_on+num_err_bins,:);
            
            all_errs(iTrial) = atan2(temp(end,2)-temp(1,2), temp(end,1) - temp(1,1));
        end
        
        dir_errs(iDir) = mean(all_errs);
    end
    clear iDir all_errs iTrial temp align_idx;
    
    % now loop along trials for each epoch
    all_errs = [];
    epoch_errs = cell(1,length(epochs));
    for iEpoch = 1:length(epochs)
        idx = find(get_trial_data_indices(td,'epoch',epochs{iEpoch}));
        [~,I] = sort([td(idx).trialId]);
        idx = idx(I); clear I;
        
        trial_errs = zeros(1,length(idx));
        for iTrial = 1:length(idx)
            dir_idx = utheta == td(idx(iTrial)).target_direction;
            
            temp = td(idx(iTrial)).vel;
            temp = temp(td(idx(iTrial)).idx_movement_on:td(idx(iTrial)).idx_movement_on+num_err_bins,:);
            
            trial_errs(iTrial) = angleDiff(atan2(temp(end,2)-temp(1,2), temp(end,1) - temp(1,1)), dir_errs(dir_idx),true,false);
        end
        
        epoch_errs{iEpoch} = trial_errs;
        all_errs = [all_errs, trial_errs];
    end
    clear iEpoch idx I trial_errs iTrial temp dir_idx;
    
    %%%%%%%%%%%%%%%%
    % Now do GPFA stuff
    num_bins = floor(1000*(align_window(2)-align_window(1))/params.bin_w)+1;
    
    % separate baseline reaches by target direction and get a mean baseline trajectory for each target
    dir_gpfa_traces = zeros(length(utheta),params.xdim,num_bins);
    for iDir = 1:length(utheta)
        idx = find(get_trial_data_indices(td,'target_direction',utheta(iDir),'epoch','BL'));
        
        all_traces = zeros(length(idx),params.xdim,num_bins);
        for iTrial = 1:length(idx)
            temp = td(idx(iTrial)).([use_array '_gpfa']);
            % gpfa might have different bin size than kinematic data
            align_idx = round(td(idx(iTrial)).(align_event)*(data_bin_w/params.bin_w));
            
            idx_1 = align_idx + floor(1000*align_window(1)/params.bin_w);
            idx_2 = align_idx + floor(1000*align_window(2)/params.bin_w);
            
            all_traces(iTrial,:,:) = temp(:,idx_1:idx_2);
        end
        
        % get mean trace
        dir_gpfa_traces(iDir,:,:) = mean(all_traces,1);
    end
    clear iDir iTrial temp align_idx all_traces;
    
    % loop through experimental epochs
    all_gpfa_diff = [];
    epoch_gpfa_diff = cell(1,length(epochs));
    for iEpoch = 1:length(epochs)
        % sort trials for each epoch by trial number (likely to already be)
        idx = find(get_trial_data_indices(td,'epoch',epochs{iEpoch}));
        [~,I] = sort([td(idx).trialId]);
        idx = idx(I); clear I;
        
        gpfa_diff = zeros(length(idx),num_bins);
        for iTrial = 1:length(idx)
            % get index for direction
            dir_idx = utheta == td(idx(iTrial)).target_direction;
            
            % get index for alignment
            align_idx = round(td(idx(iTrial)).(align_event)*(data_bin_w/params.bin_w));
            idx_1 = align_idx + floor(1000*align_window(1)/params.bin_w);
            idx_2 = align_idx + floor(1000*align_window(2)/params.bin_w);
            
            temp = td(idx(iTrial)).([use_array '_gpfa']);
            x = temp(:,idx_1:idx_2);
            
            % find euclidean distance in GPFA space for each point in window
            y = squeeze(dir_gpfa_traces(dir_idx,:,:));
            
            for t = 1:size(x,2)
                gpfa_diff(iTrial,t) = distance(x(:,t),y(:,t));
            end
        end
        
        epoch_gpfa_diff{iEpoch} = gpfa_diff;
        all_gpfa_diff = [all_gpfa_diff; gpfa_diff];
    end
    clear iEpoch idx I x temp y gpfa_diff align_idx dir_idx;
    
    results(iDate).dir_gpfa_traces = dir_gpfa_traces;
    results(iDate).epoch_gpfa_diff = epoch_gpfa_diff;
    results(iDate).all_gpfa_diff = all_gpfa_diff;
    results(iDate).dir_errs = dir_errs;
    results(iDate).epoch_errs = epoch_errs;
    results(iDate).all_errs = all_errs;
end
clear iDate dir_gpfa_traces epoch_gpfa_diff all_gpfa_diff dir_errs epoch_errs all_errs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
ymin = -0.5;
ymax = 2.2;
plotting_step_size = 10; % how many trials to average over

good_dates = find(strcmpi(sessionList(:,1),'Chewie'));

figure; hold all;
y2 = nan(length(good_dates),500);
for iDate = 1:length(good_dates)
    epoch_sizes = [150,200,min([size(results(good_dates(iDate)).epoch_gpfa_diff{3},1),150])];
    
    bl = mean(rms(results(good_dates(iDate)).epoch_gpfa_diff{1}(1:epoch_sizes(1),:),2));
    temp = [rms(results(good_dates(iDate)).epoch_gpfa_diff{1}(1:epoch_sizes(1),:),2); ...
        rms(results(good_dates(iDate)).epoch_gpfa_diff{2}(1:epoch_sizes(2),:),2); ...
        rms(results(good_dates(iDate)).epoch_gpfa_diff{3}(1:epoch_sizes(3),:),2)] - bl;
    y2(iDate,1:length(temp)) = temp;
end
clear iDate epoch_sizes bl temp;

y2 = y2./repmat(nanmean(y2,2),1,size(y2,2));
bin_vec = 1:plotting_step_size:length(y2)+1;
m = zeros(1,length(bin_vec)-1);
s = zeros(1,length(bin_vec)-1);
for iBin = 1:length(bin_vec)-1
    temp = reshape(y2(:,bin_vec(iBin):bin_vec(iBin+1)-1),1,size(y2,1)*plotting_step_size);
    m(iBin) = nanmean(temp);
    s(iBin) = nanstd(temp)./sqrt(sum(~isnan(temp)));
end
bin_vec =bin_vec(1:length(bin_vec)-1);
patch([bin_vec, fliplr(bin_vec)],[m-s, fliplr(m+s)],'b','FaceAlpha',0.5);
plot(bin_vec,m,'b-','LineWidth',2);


good_dates = find(strcmpi(sessionList(:,1),'Mihili'));

y2 = nan(length(good_dates),500);
for iDate = 1:length(good_dates)
    epoch_sizes = [150,200,min([size(results(good_dates(iDate)).epoch_gpfa_diff{3},1),150])];
    
    bl = mean(rms(results(good_dates(iDate)).epoch_gpfa_diff{1}(1:epoch_sizes(1),:),2));
    temp = [rms(results(good_dates(iDate)).epoch_gpfa_diff{1}(1:epoch_sizes(1),:),2); ...
        rms(results(good_dates(iDate)).epoch_gpfa_diff{2}(1:epoch_sizes(2),:),2); ...
        rms(results(good_dates(iDate)).epoch_gpfa_diff{3}(1:epoch_sizes(3),:),2)] - bl;
    y2(iDate,1:length(temp)) = temp;
end
clear iDate epoch_sizes bl temp;

y2 = y2./repmat(nanmean(y2,2),1,size(y2,2));
bin_vec = 1:plotting_step_size:length(y2)+1;
m = zeros(1,length(bin_vec)-1);
s = zeros(1,length(bin_vec)-1);
for iBin = 1:length(bin_vec)-1
    temp = reshape(y2(:,bin_vec(iBin):bin_vec(iBin+1)-1),1,size(y2,1)*plotting_step_size);
    m(iBin) = nanmean(temp);
    s(iBin) = nanstd(temp)./sqrt(sum(~isnan(temp)));
end
bin_vec =bin_vec(1:length(bin_vec)-1);
patch([bin_vec, fliplr(bin_vec)],[m-s, fliplr(m+s)],'r','FaceAlpha',0.5);
plot(bin_vec,m,'r-','LineWidth',2);

axis('tight');
plot([150 - rem(150,plotting_step_size),150 - rem(150,plotting_step_size)],[ymin,ymax],'k-');
plot([350 - rem(350,plotting_step_size),350 - rem(350,plotting_step_size)],[ymin,ymax],'k-');
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[ymin ymax]);
xlabel('Trial Count','FontSize',14);
ylabel('Difference','FontSize',14);
title([align_event ' ' num2str(align_window)],'FontSize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin = 0;
ymax = pi;
good_dates = find(strcmpi(sessionList(:,1),'Chewie'));

figure; hold all;
y1 = nan(length(good_dates),500);
for iDate = 1:length(good_dates)
    epoch_sizes = [150,200,min([size(results(good_dates(iDate)).epoch_errs{3},2),150])];
    
    bl = mean(results(good_dates(iDate)).epoch_errs{1}(1:epoch_sizes(1)));
    temp = angleDiff([results(good_dates(iDate)).epoch_errs{1}(1:epoch_sizes(1)), ...
        results(good_dates(iDate)).epoch_errs{2}(1:epoch_sizes(2)), ...
        results(good_dates(iDate)).epoch_errs{3}(1:epoch_sizes(3))], bl,true,true);
    y1(iDate,1:length(temp)) = temp;
end
clear iDate bl temp epoch_sizes;

y1 = y1./repmat(nanmean(y1(:,1:150),2),1,size(y1,2));
bin_vec = 1:plotting_step_size:length(y1)+1;
m = zeros(1,length(bin_vec)-1);
s = zeros(1,length(bin_vec)-1);
for iBin = 1:length(bin_vec)-1
    temp = reshape(y1(:,bin_vec(iBin):bin_vec(iBin+1)-1),1,size(y1,1)*plotting_step_size);
    m(iBin) = nanmean(temp);
    s(iBin) = nanstd(temp)./sqrt(sum(~isnan(temp)));
end
bin_vec =bin_vec(1:length(bin_vec)-1);
patch([bin_vec, fliplr(bin_vec)],[m-s, fliplr(m+s)],'b','FaceAlpha',0.5);
plot(bin_vec,m,'b-','LineWidth',2);


good_dates = find(strcmpi(sessionList(:,1),'Mihili'));
y1 = nan(length(good_dates),500);
for iDate = 1:length(good_dates)
    epoch_sizes = [150,200,min([size(results(good_dates(iDate)).epoch_errs{3},2),150])];
    
    bl = mean(results(good_dates(iDate)).epoch_errs{1}(1:epoch_sizes(1)));
    temp = angleDiff([results(good_dates(iDate)).epoch_errs{1}(1:epoch_sizes(1)), ...
        results(good_dates(iDate)).epoch_errs{2}(1:epoch_sizes(2)), ...
        results(good_dates(iDate)).epoch_errs{3}(1:epoch_sizes(3))], bl, true, true);
    y1(iDate,1:length(temp)) = temp;
end
clear iDate bl temp epoch_sizes;

y1 = y1./repmat(nanmean(y1(:,1:150),2),1,size(y1,2));
bin_vec = 1:plotting_step_size:length(y1)+1;
m = zeros(1,length(bin_vec)-1);
s = zeros(1,length(bin_vec)-1);
for iBin = 1:length(bin_vec)-1
    temp = reshape(y1(:,bin_vec(iBin):bin_vec(iBin+1)-1),1,size(y1,1)*plotting_step_size);
    m(iBin) = nanmean(temp);
    s(iBin) = nanstd(temp)./sqrt(sum(~isnan(temp)));
end
bin_vec =bin_vec(1:length(bin_vec)-1);
patch([bin_vec, fliplr(bin_vec)],[m-s, fliplr(m+s)],'r','FaceAlpha',0.5);
plot(bin_vec,m,'r-','LineWidth',2);

axis('tight');
plot([150 - rem(150,plotting_step_size),150 - rem(150,plotting_step_size)],[ymin,ymax],'k-');
plot([350 - rem(350,plotting_step_size),350 - rem(350,plotting_step_size)],[ymin,ymax],'k-');
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[ymin ymax]);
xlabel('Trial Count','FontSize',14);
ylabel('Difference','FontSize',14);


%
iEpoch = 2;
epoch_max_trials = [150,200,150];

ymin = -0.5;
ymax = 2.2;
% plotting_step_size = 20; % how many trials to average over

good_dates = find(strcmpi(sessionList(:,1),'Chewie'));

figure; hold all;
y1 = nan(length(good_dates),epoch_max_trials(iEpoch));
for iDate = 1:length(good_dates)
    epoch_sizes = [epoch_max_trials(1),epoch_max_trials(2),min([size(results(good_dates(iDate)).epoch_errs{3},2),epoch_max_trials(3)])];
    
    bl = mean(results(good_dates(iDate)).epoch_errs{1}(1:epoch_sizes(1)));
    temp = results(good_dates(iDate)).epoch_errs{iEpoch}(1:epoch_sizes(iEpoch)) - bl;
    y1(iDate,1:length(temp)) = temp;
end
clear iDate bl temp epoch_sizes;

y1 = y1./repmat(nanmean(y1,2),1,size(y1,2));
bin_vec = 1:plotting_step_size:length(y1);
m1 = zeros(size(y1,1),length(bin_vec)-1);
s1 = zeros(size(y1,1),length(bin_vec)-1);
for iBin = 1:length(bin_vec)-1
    for iDate = 1:size(y1,1)
        temp = y1(iDate,bin_vec(iBin)+1:bin_vec(iBin+1));
        m1(iDate,iBin) = nanmean(temp);
        s1(iDate,iBin) = nanstd(temp)./sqrt(sum(~isnan(temp)));
    end
end

y2 = nan(length(good_dates),epoch_max_trials(iEpoch));
for iDate = 1:length(good_dates)
    epoch_sizes = [epoch_max_trials(1),epoch_max_trials(2),min([size(results(good_dates(iDate)).epoch_errs{3},2),epoch_max_trials(3)])];
    
    bl = mean(results(good_dates(iDate)).epoch_gpfa_diff{1}(1:epoch_sizes(1)));
    temp = rms(results(good_dates(iDate)).epoch_gpfa_diff{iEpoch}(1:epoch_sizes(iEpoch),:),2) - bl;
    y2(iDate,1:length(temp)) = temp;
end
clear iDate bl temp epoch_sizes;

y2 = y2./repmat(nanmean(y2,2),1,size(y2,2));
bin_vec = 1:plotting_step_size:length(y2);
m2 = zeros(size(y2,1),length(bin_vec)-1);
s2 = zeros(size(y2,1),length(bin_vec)-1);
for iBin = 1:length(bin_vec)-1
    for iDate = 1:size(y2,1)
        temp = y2(iDate,bin_vec(iBin)+1:bin_vec(iBin+1));
        m2(iDate,iBin) = nanmean(temp);
        s2(iDate,iBin) = nanstd(temp)./sqrt(sum(~isnan(temp)));
    end
end

d1 = abs(diff(m1));
d1 = reshape(d1,size(d1,1)*size(d1,2),1);
d2 = abs(diff(m2));
d2 = reshape(d2,size(d2,1)*size(d2,2),1);
plot(d1,d2,'bd','LineWidth',2);

% fit line and plot it
[b,~,~,~,s] = regress(d2,[ones(size(d1)) d1]);
plot(d1,b(1)+b(2)*d1,'b-','LineWidth',2);

good_dates = find(strcmpi(sessionList(:,1),'Mihili'));

y1 = nan(length(good_dates),epoch_max_trials(iEpoch));
for iDate = 1:length(good_dates)
    epoch_sizes = [epoch_max_trials(1),epoch_max_trials(2),min([size(results(good_dates(iDate)).epoch_errs{3},2),epoch_max_trials(3)])];
    
    bl = mean(results(good_dates(iDate)).epoch_errs{1}(1:epoch_sizes(1)));
    temp = results(good_dates(iDate)).epoch_errs{iEpoch}(1:epoch_sizes(iEpoch)) - bl;
    y1(iDate,1:length(temp)) = temp;
end
clear iDate bl temp epoch_sizes;

y1 = y1./repmat(nanmean(y1,2),1,size(y1,2));
bin_vec = 1:plotting_step_size:length(y1);
m1 = zeros(size(y1,1),length(bin_vec)-1);
s1 = zeros(size(y1,1),length(bin_vec)-1);
for iBin = 1:length(bin_vec)-1
    for iDate = 1:size(y1,1)
        temp = y1(iDate,bin_vec(iBin)+1:bin_vec(iBin+1));
        m1(iDate,iBin) = nanmean(temp);
        s1(iDate,iBin) = nanstd(temp)./sqrt(sum(~isnan(temp)));
    end
end

y2 = nan(length(good_dates),epoch_max_trials(iEpoch));
for iDate = 1:length(good_dates)
    epoch_sizes = [epoch_max_trials(1),epoch_max_trials(2),min([size(results(good_dates(iDate)).epoch_errs{3},2),epoch_max_trials(3)])];
    
    bl = mean(results(good_dates(iDate)).epoch_gpfa_diff{1}(1:epoch_sizes(1)));
    temp = rms(results(good_dates(iDate)).epoch_gpfa_diff{iEpoch}(1:epoch_sizes(iEpoch),:),2) - bl;
    y2(iDate,1:length(temp)) = temp;
end
clear iDate bl temp epoch_sizes;

y2 = y2./repmat(nanmean(y2,2),1,size(y2,2));
bin_vec = 1:plotting_step_size:length(y2);
m2 = zeros(size(y2,1),length(bin_vec)-1);
s2 = zeros(size(y2,1),length(bin_vec)-1);
for iBin = 1:length(bin_vec)-1
    for iDate = 1:size(y2,1)
        temp = y2(iDate,bin_vec(iBin)+1:bin_vec(iBin+1));
        m2(iDate,iBin) = nanmean(temp);
        s2(iDate,iBin) = nanstd(temp)./sqrt(sum(~isnan(temp)));
    end
end

d1 = abs(diff(m1));
d1 = reshape(d1,size(d1,1)*size(d1,2),1);
d2 = abs(diff(m2));
d2 = reshape(d2,size(d2,1)*size(d2,2),1);
plot(d1,d2,'ro','LineWidth',2);

[b,~,~,~,s] = regress(d2,[ones(size(d1)) d1]);
plot(d1,b(1)+b(2)*d1,'r-','LineWidth',2);



%%
for iEpoch = 1:length(epochs)
    [y,g] = deal([]);
    for iDate = 1:length(good_dates)
        switch iEpoch
            case 1
                tempy = y1(iDate,1:150);
                tempg = 1:150;
            case 2
                tempy = y1(iDate,151:350);
                tempg = 1:200;
            case 3
                tempy = y1(iDate,351:500);
                tempg = 1:150;
        end
        
        idx = isnan(tempy);
        tempy = tempy(~idx);
        tempg = tempg(~idx);
        
        y = [y, tempy];
        g = [g,tempg];
    end
    
    p(iEpoch) = anovan(y',g','continuous',1,'display','off');
end
clear iEpoch iDate tempy tempg y g;

p

%%
good_dates = [7,8,9,10,11,12,14,17];
good_dates = [13,15,16,18,19,20];
