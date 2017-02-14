%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = parseFileByTrial(data, params)
%
%   Wrapper function to create trial_data structs. Currently only supports
% CDS.
%
% INPUTS:
%   data     : a CDS file or BDF file
%   params   : a struct containing parameters
%     .meta         : a struct with a field for each meta parameter you want attached
%                       to this file. This can handle any arbitrary information!
%     .trialResults : which reward codes to use ('R','A','F','I')
%     .binSize      : default 0.01 sec
%     .extraTime    : [time before, time after] beginning and end of trial (default [0.5 0.3] sec)
%
% OUTPUTS:
%   trial_data : the struct! Huzzah!
%
% TIPS:
%   I recommend the following extra function calls afterwards
%       trial_data = getMoveOnsetAndPeak(trial_data);
%       trial_data = pruneBadTrials(trial_data);
%       trial_data = getCommonUnits(trial_data);
%
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = parseFileByTrial(data,params,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1, params = []; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(data,'commonDataStructure')
    trial_data = parseFileByTrial_cds(data,params);
else
    error('BDF not currently supported.');
    % when someone is ready to implement this, can write the following func
    %   trial_data = parseFileByTrial_bdf(data,params);
end

% ensure standard field order
trial_data = reorderTDfields(trial_data);
