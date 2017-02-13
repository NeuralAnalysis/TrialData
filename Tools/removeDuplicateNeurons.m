function [trial_data,perc_coinc] = removeDuplicateNeurons(trial_data)
% function to check for shunts or duplicate neurons

all_spikes = cat(1,trial_data.M1_spikes);

perc_coinc = NaN(size(all_spikes,2));
for unit1 = 1:size(all_spikes,2)
    for unit2 = unit1:size(all_spikes,2)
        perc_coinc(unit1,unit2) = 100*sum(all_spikes(:,unit1) & all_spikes(:,unit2))/size(all_spikes,1);
    end
end

all_dist = reshape(perc_coinc,numel(perc_coinc),1);

% use all_dist to identify cells that have outlier coincidence and remove

% IDEAS FOR EXCLUSION CRITERIA
%   1) build "master" distribution using TONS of sessions and save it here
%   permanently (fit a distribution to it?)
%   2) use all pairs in a session to look for outlier pairs
%   3) give an arbitrary cutoff

root_dir = '/Users/mattperich/Data/TrialDataFiles/';
fn=dir([root_dir '*.mat']);
all_dist = [];
for i = 1:length(fn)
    load([root_dir fn(i).name]);
    
    all_spikes = cat(1,trial_data.M1_spikes);

perc_coinc = NaN(size(all_spikes,2));
for unit1 = 1:size(all_spikes,2)
    for unit2 = unit1:size(all_spikes,2)
        perc_coinc(unit1,unit2) = 100*sum( (all_spikes(:,unit1) > 0) & ...
                                            (all_spikes(:,unit2) > 0) & ...
                                            (all_spikes(:,unit1) == all_spikes(:,unit2))) ...
                                            /size(all_spikes,1);
    end
end
all_dist = [all_dist; reshape(perc_coinc,numel(perc_coinc),1)];
end

dist_prctile = zeros(1,100);
for i = 1:100
    
end


