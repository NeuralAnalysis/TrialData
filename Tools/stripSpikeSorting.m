function trial_data = stripSpikeSorting(trial_data)
% as the name implies, this strips away the spike sorting

fn_array = getTDfields(trial_data,'arrays');

for iArray = 1:length(fn_array)
    sg = trial_data(1).([fn_array{iArray} '_unit_guide']);
    elecs = unique(sg(:,1));
    
    % initialize the temporary field
    for trial = 1:length(trial_data)
        trial_data(trial).temp_spikes = NaN(size(trial_data(trial).([fn_array{iArray} '_spikes']),1),length(elecs));
    end
    
    for e = 1:length(elecs)
        idx = sg(:,1) == elecs(e);
        
        temp = trial_data(trial).([fn_array{iArray} '_spikes']);
        trial_data(trial).temp_spikes(:,e) = sum(temp(:,idx),2);
    end
    
    
    % lastly, make and distribute new spike guide and finalize fields
    sg_new = [elecs, zeros(length(elecs),1)];
    for trial = 1:length(trial_data)
        trial_data(trial).([fn_array{iArray} '_spikes']) = trial_data(trial).temp_spikes;
        trial_data(trial).([fn_array{iArray} '_unit_guide']) = sg_new;
    end
    trial_data = rmfield(trial_data,'temp_spikes');
end
