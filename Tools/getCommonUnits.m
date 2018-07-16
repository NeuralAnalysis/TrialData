%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = getCommonUnits(trial_data)
%
%   Looks across all trials of a trial_data struct and ensure the units
% are all common. This is mostly a sanity check/housekeeping thing.
%
% INPUTS:
%   trial_data : the struct
%
% OUTPUTS:
%   trial_data : same struct, but with unmatched units removed
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = getCommonUnits(trial_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end

% figure out which arrays are here
arrays = getTDfields(trial_data,'arrays');

for idx_array = 1:length(arrays)
    array = arrays{idx_array};
    
    master_sg = trial_data(1).([array '_unit_guide']);
    % get global unit guide
    for i = 2:length(trial_data)
        bad_idx = checkUnitGuides(master_sg,trial_data(i).([array '_unit_guide']));
        if ~isempty(bad_idx)
%             disp('Found a bad unit...');
            master_sg = setdiff(master_sg, bad_idx, 'rows');
        end
    end
    
    % now trim out missing units
    for i = 1:length(trial_data)
        [~, idx] = checkUnitGuides(trial_data(i).([array '_unit_guide']),master_sg);
        temp = trial_data(i).([array '_spikes']);
        trial_data(i).([array '_spikes']) = temp(:,~idx);
        temp = trial_data(i).([array '_unit_guide']);
        trial_data(i).([array '_unit_guide']) = temp(~idx,:);
    end
end

