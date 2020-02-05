function out = aliasTDfields(trial_data,params)
% ALIASTDFIELDS copies and aliases fields of trial_data accoring to params
% Inputs:
%   trial_data: the trial_data struct
%   params: params struct
%       .alias_list - cell array of names, where first column is field to be renamed
%           and second column is the new name
%       .overwrite_fields - boolean value for whether old fields should be overwritten if they exist (default: false)
%       .remove_old_fields - boolean value for whether to remove old fields (default: true)
%
% Written by Raeed Chowdhury, Feb 2020

alias_list = {};
overwrite_fields = false;
remove_old_fields = true;
assignParams(who,params)

out = trial_data;

assert(size(alias_list,2)==2 || isempty(alias_list),'Alias list must be a two column list of fields to rename and their new names')
for i = 1:size(alias_list,1)
    if isfield(trial_data,alias_list{i,1})
        if overwrite_fields || ~isfield(out,alias_list{i,2})
            [out.(alias_list{i,2})] = trial_data.(alias_list{i,1});
            if remove_old_fields
                out = rmfield(out,alias_list{i,1});
            end
        else
            warning('Field named %s exists, but overwrite flag is not set. Field was not overwritten',alias_list{i,2})
        end
    else
        warning('Field named %s does not exist',alias_list{i,1})
    end
end
out = reorderTDfields(out);
