function cosDelta = cosReflectionFromPointOnSphere(lon, lat, R, d, phase_angle)
% cosReflectionFromPointOnSphere
% Returns the cosine of the local reflection angle between the surface normal
% at a given point (lon, lat) on a sphere of radius R and the vector toward
% an observer located at [d*cos(phase_angle), d*sin(phase_angle), 0].
%
% Inputs:
%   lon         - Longitude φ (radians)
%   lat         - Latitude λ (radians)
%   R           - Radius of the sphere
%   d           - Observer distance from center of sphere (in equatorial plane)
%   phase_angle - Phase angle α (radians)
%
% Output:
%   cosDelta    - Cosine of the local reflection angle (dimensionless)
%
% --------  DERIVATION (updated)  ------------------------------------------
% • Surface point p (also surface normal n = p/R):
%       p = R [ cosλ cosφ,  cosλ sinφ,  sinλ ]                        (1)
%
% • Observer position o:
%       o = [ d cosα,  d sinα,  0 ]                                   (2)
%
% • Vector from surface point to observer:
%       v = o − p                                                    (3)
%
% • Surface normal n = p / R                                         (4)
%
% • cosΔ ≡ cos(θ_r) = (n · v) / |v|                                  (5)
%
% • Dot product n·v using (1)–(4):
%       n·v = (1/R) * (o · p − |p|²)                                 (6)
%            = d cosλ cos(φ−α) − R                                   (7)
%
% • Magnitude of vector v:
%       |v|² = |o|² + |p|² − 2 o·p                                   (8)
%            = d² + R² − 2 d R cosλ cos(φ−α)                         (9)
%
% • Final formula:
%       cosΔ = [d cosλ cos(φ−α) − R] / sqrt(d² + R² − 2 d R cosλ cos(φ−α))   (10)
% --------------------------------------------------------------------------

d_point = distanceToPointOnSphere(lon, lat, R, d, phase_angle);

cosDelta = (d.*cos(lat).*cos(lon - phase_angle) - R) ./ d_point;

end
