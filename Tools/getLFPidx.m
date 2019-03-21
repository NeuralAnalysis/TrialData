function [idx,freqs,chans] = getLFPidx(td,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [idx,freqs,chans] = getLFPidx(td, chans, freqs)
%
%   This function will query the lfp_guide and return indices.
% 
% INPUTS:
%   td               :  the struct
%   params           :  params struct
%       .array       :  which array to use
%       .channels    :  which channels to return. Use [] or 'all' for all
%       .freq_bands  :  which frequency ranges [LOW, HIGH]
%                           use [] or 'all' for all
%
%   OR
%
%   Instead of a params struct, ...'NAME',VALUE... pairs
%
% OUTPUTS:
%   idx       :  
%   freqs     :  vector listing all frequency bands
%   chans     :  vector listing all channels in the  lfp
%
% Written by Matt Perich.  Updated March 2019.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array       =  '';
channels    =  'all';
freq_bands  =  'all';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    if nargin == 2 % it's a params struct
        params = varargin{1};
        if isstruct(params)
            assignParams(who,params);
        else
            error('Could not parse input to getLFPidx. Should either be a params struct or ...''NAME'',VALUE... pairs.');
        end
    else % it's name value pairs
        if mod(nargin,2) ~= 0 % make sure there are pairs
            error('Could not parse input to getLFPidx. Should either be a params struct or ...''NAME'',VALUE... pairs.');
        else
            params = struct();
            for iVar = 1:2:length(varargin)
                params.(varargin{iVar})  =  varargin{iVar+1};
            end
            assignParams(who,params);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(array)
    % first check if there is more than one LFP
    fn = getTDfields(td,'lfp');
    if isempty(fn) % no LFP, no need for this function!
        disp('ERROR: getTDfields: No LFP field found!');
        idx = [];
        freqs = [];
        chans = [];
        return;
    elseif length(fn) > 1
        disp('ERROR: getTDfields: Multiple  LFP fields found...  you need to specify which array.');
        idx = [];
        freqs = [];
        chans = [];
        return;
    elseif length(fn) == 1
        array = getTDfields(td,'arrays');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lfp_guide = td(1).([array '_lfp_guide']);
freqs = unique(lfp_guide(:,2:3),'rows');
chans = unique(lfp_guide(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(channels) || (ischar(channels) && strcmpi(channels,'all'))
    channels = chans;
end

if isempty(freq_bands) || (ischar(freq_bands) && strcmpi(freq_bands,'all'))
    freq_bands = freqs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = find(ismember(lfp_guide(:,1),channels) & ...
    ismember(lfp_guide(:,2:3),freq_bands,'rows'));


