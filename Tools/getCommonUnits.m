function trial_data = getCommonUnits(trial_data)
% will look across all trials of a trial_data struct and ensure the units
% are all common

% figure out which arrays are here
fn = fieldnames(trial_data);
fn = fn(cellfun(@(x) ~isempty(x),strfind(fieldnames(trial_data),'_spikes')));
arrays = strrep(fn,'_spikes','')';

for idx_array = 1:length(arrays)
    array = arrays{idx_array};
    
    master_sg = trial_data(1).([array '_unit_guide']);
    % get global unit guide
    for i = 2:length(trial_data)
        bad_idx = checkUnitGuides(master_sg,trial_data(i).([array '_unit_guide']));
        if ~isempty(bad_idx)
            disp('Found a bad unit...');
            master_sg = setdiff(master_sg, bad_idx, 'rows');
        end
    end
    
    % get rid of things that aren't utah recordings
    master_sg( master_sg(:,1) > 96, :) = [];
    
    % now trim out missing units
    for i = 1:length(trial_data)
        [~, idx] = checkUnitGuides(trial_data(i).([array '_unit_guide']),master_sg);
        temp = trial_data(i).([array '_spikes']);
        trial_data(i).([array '_spikes']) = temp(:,~idx);
        temp = trial_data(i).([array '_unit_guide']);
        trial_data(i).([array '_unit_guide']) = temp(~idx,:);
    end
end

