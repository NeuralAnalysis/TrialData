function [idx,freqs,chans] = getLFPidx(td,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [idx,freqs,chans] = getLFPidx(td, chans, freqs)
%
%   This function 
% 
% INPUTS:
%   td               :  the struct
%   params           :  params struct
%       .array       :  which array to use
%       .channels    :  which channels to return. Use [] or 'all' for all
%       .freq_bands  :  which frequency ranges [LOW, HIGH]
%                           use [] or 'all' for all
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
    assignParams(who,params);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(array)
    % first check if there is more than one LFP
    fn = getTDfields(td,'lfp');
    if isempty(fn)
        error('No LFP field found');
    elseif length(fn) > 1
        error('Multiple  LFP fields found...  you need to specify which array.');
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


