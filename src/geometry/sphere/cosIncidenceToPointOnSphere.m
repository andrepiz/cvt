function cosDelta = cosIncidenceToPointOnSphere(lon, lat, R, d)
% cosIncidenceToPointOnSphere
% Returns the cosine of the local incidence angle between the surface normal
% at a given point (lon, lat) on a sphere of radius R and the vector toward
% a light located at [d, 0, 0].
%
% Inputs:
%   lon         - Longitude φ (radians)
%   lat         - Latitude λ (radians)
%   R           - Radius of the sphere
%   d           - Light distance from center of sphere (in equatorial plane)
%
% Output:
%   cosDelta    - Cosine of the local incidence angle (dimensionless)


cosDelta = cosReflectionFromPointOnSphere(lon, lat, R, d, 0);

end
