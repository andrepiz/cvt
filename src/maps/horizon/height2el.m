function [elA, elB] = height2el(heightA, lonA, latA, heightB, lonB, latB)
% height2el
% Computes the elevation angle (in radians) from a set of points A located at
% (lonA, latA) with height heightA to a set of points B (lonB, latB)
% with heights given by heightB.
%
% Inputs:
%   heightA - Height(s) of observer point(s) A above the reference sphere (same units as heightB)
%   lonA    - Longitude(s) of observer point(s) A (in radians)
%   latA    - Latitude(s) of observer point(s) A (in radians)
%   heightB - Height(s) of target point(s) B above the reference sphere
%   lonB    - Longitude(s) of target point(s) B (in radians)
%   latB    - Latitude(s) of target point(s) B (in radians)
%
% Outputs:
%   elA - Elevation angle(s) from point A to point B (in radians)
%   elB - Elevation angle(s) from point B to point A (in radians)
%
% Notes:
% - Positive elevation angles indicate the target is above the local horizon.
% - Negative values indicate the target is below the horizon.

% Compute great-circle angular distance gamma between points A and B
cosGamma = sin(latB) .* sin(latA) + cos(latB) .* cos(latA) .* cos(lonB - lonA);

% Calculate elevation angle using spherical triangle relationships
coElFac = sin(acos(cosGamma)) ./ sqrt(heightA.^2 + heightB.^2 - 2*heightA.*heightB.*cosGamma);
sinCoElA = heightB .* coElFac;
elA = pi/2 - asin(max(-1, min(1, sinCoElA)));

% Reverse the angle of those points that are located below the tangent plane
ix_elA_reverse = heightB .* cosGamma < heightA;
elA(ix_elA_reverse) = -elA(ix_elA_reverse);

if nargout > 1
    sinCoElB = heightA .* coElFac;
    elB = pi/2 - asin(max(-1, min(1, sinCoElB)));
    ix_elB_reverse = heightA .* cosGamma < heightB;
    elB(ix_elB_reverse) = -elB(ix_elB_reverse);
end

end
