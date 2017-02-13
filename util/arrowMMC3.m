% useage:
%   h = arrowMMC(prevPoint, point, nextPoint)
%       where each of the three points is [x,y].  The three points are used to set the orientation
%       of the arrow so that it looks good when plotted on top of a curve, and roughly follows the
%       tangent of the curvature.
%
%       If you don't want to use the third point (e.g., if there isn't one) then use:
%   h = arrowMMC(prevPoint, point, [])
%
%       You can also specify these additional arguments (must be done in order, but you don't have to use all)
%   h = arrowMMC(prevPoint, point, nextPoint, sizeArrow, axisRange, faceColor, edgeColor)
%       Most importantly, axisRange tells arrowMMC how big to make the arrow, and allows it to have
%       appropriately scaled proportions.  If you don't supply this, arrowMMC will get it using
%       'gca'.  But, if the axes are then later rescaled (e.g. due to more points plotting) things
%       will get wonky.  So either supply 'axisRange', or make sure the axes don't change after
%       plotting the arrow.
%
%       A reasonable starting size for the arrow is 6.
%
function h = arrowMMC3(prevPoint, point, nextPoint, sizeArrow, axisRange, faceColor, edgeColor)


roughScale = 0.004;  % setting this empirically so that 'sizeArrow' works roughly like 'markerSize'
xVals = [0 -1.5 4.5 -1.5 0];
yVals = [0 2 0 -2 0];
zVals = [0 2 0 -2 0];


%% do some parsing of the inputs, and use defaults if not provided
if isempty(nextPoint)
    nextPoint = point + point-prevPoint;
end

if nargin<4, sizeArrow = 6; end

if nargin<5
    xRange = range(get(gca,'XLim'));
    yRange = range(get(gca,'YLim'));
    zRange = range(get(gca,'ZLim'));
else
    xRange = range(axisRange(1:2));
    yRange = range(axisRange(3:4));
    zRange = range(axisRange(5:6));
end

if nargin<6, faceColor = [0 0 0]; end
if nargin<7, edgeColor = faceColor; end

%% do a bit of scaling
mxX = max(xVals);
xVals = roughScale*sizeArrow * xVals/mxX * xRange;
yVals = roughScale*sizeArrow * yVals/mxX * yRange;
zVals = roughScale*sizeArrow * zVals/mxX * zRange;

%% now do the rotation

vector = nextPoint - prevPoint;
theta = atan2(vector(2),vector(1));
% rotM = [cos(theta) -sin(theta); sin(theta), cos(theta)];
rotM = [cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];

newVals = rotM*[xVals; yVals; zVals];
xVals = newVals(1,:);
yVals = newVals(2,:);
zVals = newVals(3,:);

%% now plot
xVals = xVals + point(1);
yVals = yVals + point(2);
zVals = zVals + point(3);
h = fill3(xVals, yVals,zVals, faceColor);
set(h, 'edgeColor', edgeColor');

