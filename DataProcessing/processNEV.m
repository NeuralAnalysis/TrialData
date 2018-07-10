

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = processNEV(filename,signal_info)
% process .nev files
%
% some parameters
data_type      = 'spikes'; % 'spikes', 'comments' for now
save_NEV_mat   = false;
reload_NEV     = false;
read_waveforms = false; % (not implemented) whether to include waveforms for spikes. Big file size increase
spiking_chans  = 1:96;
exclude_units  = 255; % sort id of units to exclude
strip_sort     = false;  % whether to ignore unit sort codes
assignParams(who,signal_info.params); % overwrite parameters

% here is where I can prepare the inputs if I want
%   but the openNEV has a super-shit design and you can't have variable
%   inputs without building a string eval command. Dumbasses.
[~,f,~] = fileparts(filename);
if exist([f,'.mat'],'file') && ~reload_NEV
    load([f,'.mat'],'NEV');
else
    if save_NEV_mat
        NEV = openNEV_td(filename);
    else
        NEV = openNEV_td(filename,'nosave');
    end
end

if isempty(NEV)
    error('NEV not found.');
end
%%%%%%%%%%%
% TO DO
%   Add support for captured digital events
%%%%%%%%%%%


%%%%%%
% parse NEV format into something easier to use to get spiking data into TD
count = 0;
out = struct( ...
    't', 0:1/double(NEV.MetaTags.SampleRes):NEV.MetaTags.DataDurationSec, ...
    'labels',zeros(length(spiking_chans),2), ...
    'data',{{}}, ...
    'wf',{{}}); % for future waveform implementation

switch lower(data_type)
    
    case 'spikes'
        
        for iElec = 1:length(spiking_chans)
            chan_idx = NEV.Data.Spikes.Electrode == spiking_chans(iElec);
            
            % if no spikes are found, we just need to make sure it gets
            % filled with something
            if strip_sort || sum(chan_idx) == 0
                found_units = 0;
            else
                found_units = unique(NEV.Data.Spikes.Unit(chan_idx));
                % exclude based on unit id, if sorted
                found_units = setdiff(found_units,exclude_units);
            end
            
            for iUnit = 1:length(found_units)
                count = count + 1;
                if strip_sort
                    unit_idx = chan_idx;
                else
                    unit_idx = chan_idx & NEV.Data.Spikes.Unit == found_units(iUnit);
                end
                
                out.labels(count,:) = [spiking_chans(iElec), found_units(iUnit)];
                out.data{count} = double(NEV.Data.Spikes.TimeStamp(unit_idx))/NEV.MetaTags.TimeRes;
                if read_waveforms
                    out.wf{count} = NEV.Data.Spikes.Waveform(:,unit_idx);
                end
                
            end
        end
        
    case 'comments'
        % get list of text
        nev_text = reshape( ...
            NEV.Data.Comments.Text, ...
            size(NEV.Data.Comments.TimeStamp,2), ...
            numel(NEV.Data.Comments.Text)/size(NEV.Data.Comments.TimeStamp,2));
        
        labels = unique(nev_text,'rows');
        labels = mat2cell(labels,ones(1,size(labels,1)),size(labels,2));
        
        data = cell(1,length(labels));
        for iText = 1:length(labels)
            idx = startsWith(nev_text,labels{iText});
            data{iText} = NEV.Data.Comments.TimeStampSec(idx);
        end
        
        % now process the labels to make them more suited to names. NEV
        % pads the labels with blank spaces
        for iText = 1:length(labels)
            temp = labels{iText};
            temp = temp(int16(temp) > 0);
            labels{iText} = temp;
        end
        
        out.labels = labels;
        out.data = data;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
