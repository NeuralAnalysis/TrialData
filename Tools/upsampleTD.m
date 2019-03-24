function trial_data = upsampleTD(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = upsampleTD(trial_data, params)
%
%   Upsamples TD via interp.
%
% INPUTS:
%   trial_data   :  the struct
%   params       :  params struct
%       .upsample_factor  :  integer factor to upsample
%           note can just pass  in the integer without a params struct if
%           you want.
%
% OUTPUTS:
%   trial_data   :   the upsampled struct
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
upsample_factor = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if isstruct(params)
        assignParams(who,params);
    else % default to
        upsample_factor = params;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn_time = getTDfields(trial_data,'time');
fn_idx  = getTDfields(trial_data,'idx');

new_bin_size = trial_data(1).bin_size /  upsample_factor;

for trial = 1:length(trial_data)
    % first upsample time signals
    for iSig = 1:length(fn_time)
        temp = trial_data(trial).(fn_time{iSig});
        temp_us = zeros(size(temp,1)*upsample_factor,size(temp,2));
        for iCol = 1:size(temp,2)
            temp_us(:,iCol) = interp(temp(:,iCol),upsample_factor);
        end
        trial_data(trial).(fn_time{iSig}) = temp_us;
    end
    
    % now upsample idx_ fields
    for iSig = 1:length(fn_idx)
        trial_data(trial).(fn_idx{iSig}) = trial_data(trial).(fn_idx{iSig}) * upsample_factor;
    end
    
    % add the new bin size
    trial_data(trial).bin_size = new_bin_size;
end



