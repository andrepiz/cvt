function omega = solidAngleToPointOnSphere(lon, lat, R, d, phase_angle, A)
% SOLIDANGLETOPOINTONSPHERE  Find the solid angle between a point at longitude 
% and latitude coordinates on a sphere and an observer placed at a given 
% phase angle with given observing area.
%
%   omega = solidAngleToPointOnSphere(lon, lat, R, d, phase_angle, A)
%
% Inputs
% ------
%   lon          Longitude of surface point S (radians, range: 0 to 2π).
%   lat          Latitude of surface point S  (radians, range: –π/2 to π/2).
%   R            Radius of the sphere (positive scalar or array).
%   d            Distance from origin to observer point P (positive scalar).
%   phase_angle  Azimuthal angle of observer P in XY-plane (radians, 0 to 2π).
%   A            Collecting (pupil) area of the observer (in same units as R²).
%
% Geometry
% --------
%   • Sphere is centered at origin O = (0, 0, 0).
%   • Observer point P lies in the XY-plane:
%         P = (d * cos(phase_angle), d * sin(phase_angle), 0)
%   • Surface point S is located on the sphere:
%         S = (R * cos(lat) * cos(lon), 
%              R * cos(lat) * sin(lon), 
%              R * sin(lat))
%
% Solid Angle Calculation
% -----------------------
%   The solid angle ω subtended by the observer’s aperture A, as seen from
%   surface point S, is:
%       ω = A / |P - S|²
%   assuming the collecting area A is perpendicular to the line of sight
%   from S to P.
%
% Author:  Andrea Pizzetti
% Updated: 25-Jun-2025

cos_delta = cosBearingToPointOnSphere(lon, lat, R, d, phase_angle);
d_point = distanceToPointOnSphere(lon, lat, R, d, phase_angle);

omega = A * cos_delta./ d_point.^2;

end
