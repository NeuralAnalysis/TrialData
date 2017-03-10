%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function theta = binAngles(theta,angle_bin_size)
%
%   Bins angles into equal increments of angle_bin_size.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta = bin_angles(theta,angle_bin_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = minusPi2Pi(theta);
theta = round(theta./angle_bin_size).*angle_bin_size;

% Now do some checks to make sure the values make sense
% -pi and pi are the same thing
if length(unique(theta)) > int16(2*pi/angle_bin_size)
    % probably true that -pi and pi both exist
    utheta = unique(theta);
    if utheta(1)==-utheta(end)
        % almost definitely true that -pi and pi both exist
        theta(theta==utheta(1)) = utheta(end);
    elseif abs(utheta(1)) > abs(utheta(end))
        theta(theta==utheta(1)) = -utheta(1);
        % probably means that -pi instead of pi
    else
        disp('Something fishy is going on with this binning...')
    end
end