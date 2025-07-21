function elevation = elevationFromPointOnSphere(lon, lat, R, d, phase_angle)
% ELEVATIONFROMPOINTONSPHERE
% Returns the elevation angle (in radians) between the local tangent plane
% at a point (lon, lat) on a sphere and the vector toward an observer located at
% [d*cos(phase_angle), d*sin(phase_angle), 1].
%
% Inputs:
%   lon         - Longitude φ (radians)
%   lat         - Latitude λ (radians)
%   R           - Radius of the sphere
%   d           - Observer distance from center of sphere (in equatorial plane)
%   phase_angle - Phase angle α (radians)
%
% Output:
%   elevation   - Elevation angle from local tangent plane (radians)

% Convert to elevation angle (sin(θₑ) = cos(θᵣ))
elevation = asin(cosReflectionFromPointOnSphere(lon, lat, R, d, phase_angle));

end
