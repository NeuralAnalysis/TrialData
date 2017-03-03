function adiff = angleDiff(angle1,angle2,useRad,preserveSign)
% ANGLEDIFF finds absolute value of difference in angle
% 
% INPUTS:
%   angle1: an angle or vector of angles (DOESN'T WORK WITH MATRICES YET)
%   angle2: another angle or vector of angles. Same size as angle1
%   useRad: (boolean) if true, inputs are in radians
%   preserveSign: (boolean) if true, will not return absolute value
% 
% OUTPUTS:
%   adiff: element-wise difference between the angles
% 
%%%%%%
% written by Matt Perich; last updated April 2014
%%%%%%

if nargin < 4
    preserveSign = true;
    if nargin < 3 %default to degrees
        useRad = true;
    end
end

if useRad
    a = pi;
else
    a = 180;
end

% Find the difference
adiff = angle2 - angle1;

% If greater than 180, subtract from 360 to get magnitude
for i = 1:length(adiff)
    % if it's more than 360, subtract down
    if abs(adiff(i)) > 2*a
        adiff(i) = mod(adiff(i),sign(adiff(i))*2*a);
    end
    
    if adiff(i) > a
        adiff(i) = adiff(i) - 2*a;
    elseif adiff(i) < -a
        adiff(i) = adiff(i) + 2*a;
    end
end

% preserve the sign... more counterclockwise is positive
if ~preserveSign
    adiff = abs(adiff);
end