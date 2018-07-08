function [badUnits, varargout] = checkUnitGuides(varargin)
% CHECKUNITS Checks to ensure units are consistent across all inputs. Unit
% IDs are expected to be defined by the first column.
%
% each varargin is a spike guide, but has option to pass in a cell array of
% spike guides as well. So, if varargin is a single cell, assumes it is a
% cell array of guides.
%
% Example:
%   sg1 = [1 1; 1 2; 3 1; 4 1; ......]; sg2 = [ .......
%   badUnits = checkUnitGuides(sg1,sg2,sg3);
%       or
%   allSG = {sg1, sg2, sg3};
%   badUnits = checkUnitGUides( allSG );
%
% From this, to find the set that is common to all of them
%   sg = setdiff(sg1, badUnits, 'rows');

if iscell(varargin{1})
    varargin = varargin{1};
end

badUnits = [];
for i = 2:length(varargin)

    temp1 = varargin{i-1};
    temp2 = varargin{i}; 

    % If either is empty they can't really be compared...
    if ~isempty(temp1) && ~isempty(temp2)
        badUnits = [badUnits; setxor(temp1,temp2,'rows')];
    end
end

badUnits = unique(badUnits,'rows');

% return indices of rows with badUnits for each guide
for i = 1:length(varargin)
    varargout{i} = ismember(varargin{i},badUnits,'rows');
end