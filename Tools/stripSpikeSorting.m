%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = stripSpikeSorting(trial_data)
%
%   trial_data = stripSpikeSorting(trial_data)
%   This function, as the name implies, strips away the sort codes for all
% of the spiking arrays. Essentially, summing all of the spike counts in
% each bin for all units identified on a single electrode.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = stripSpikeSorting(trial_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(trial_data), error('First input must be trial_data struct!'); end
% get the arrays that are present
fn_array = getTDfields(trial_data,'arrays');

% loop along the arrays
for iArray = 1:length(fn_array)
    sg = trial_data(1).([fn_array{iArray} '_unit_guide']);
    elecs = unique(sg(:,1));
    
    % initialize the temporary field
    for trial = 1:length(trial_data)
        trial_data(trial).temp_spikes = NaN(size(trial_data(trial).([fn_array{iArray} '_spikes']),1),length(elecs));
    end
    
    % sum the spikes for all units on each electrode
    for e = 1:length(elecs)
        idx = sg(:,1) == elecs(e);
        
        for trial = 1:length(trial_data)
            temp = trial_data(trial).([fn_array{iArray} '_spikes']);
            % store result in the temporary field
            trial_data(trial).temp_spikes(:,e) = sum(temp(:,idx),2);
        end
    end
    
    
    % lastly, make and distribute new spike guide and finalize fields
    sg_new = [elecs, zeros(length(elecs),1)];
    for trial = 1:length(trial_data)
        trial_data(trial).([fn_array{iArray} '_spikes']) = trial_data(trial).temp_spikes;
        trial_data(trial).([fn_array{iArray} '_unit_guide']) = sg_new;
    end
    % get rid of this useless old thing
    trial_data = rmfield(trial_data,'temp_spikes');
end
