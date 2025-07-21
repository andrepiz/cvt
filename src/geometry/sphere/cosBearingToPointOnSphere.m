function cosBearing = cosBearingToPointOnSphere(lon, lat, R, d, phase_angle)
% ANGLESPHEREPOINTPHASE  Find the cosine of the angle delta between the vector linking an observer to 
% the center of a sphere of radius R placed at distance d and the vector
% linking the observer with a point at longitude and latitude coordinates
% on the same sphere. The observer is placed at a given phase angle from
% the sphere. 
%
%   cosBearing = cosBearingToPointOnSphere(lon, lat, R, d, phase_angle) returns the
%   cosine of the angle delta between
%        • v2 = S – P  (from P to a point S on the sphere surface)
%        • v1 = C – P  (from P to the origin C)
%
%   Inputs
%   ------
%   lon          Longitude of S (radians,  0 … 2π).
%   lat          Latitude of S  (radians, –π/2 … π/2).
%   R            Sphere radius (positive scalar or array).
%   d            Distance |OP| from the origin to point P (positive).
%   phase_angle  Azimuth of P in the XY-plane, measured from +X (radians, 0 … 2π).
%
% Geometry
% --------
%   • Sphere centre C is at the origin O = (0,0,0).
%   • Point P sits in the XY-plane:
%         P = (d cos(phase_angle), d sin(phase_angle), 0).
%   • Surface point S in standard spherical coordinates:
%         S = (R cos(lat) cos(lon), R cos(lat) sin(lon), R sin(lat)).
%
% Derivation
% ------------------
%   v1 = –P,  v2 = S – P.
%   v1·v2 = d² – dR cos(lat) cos(lon – phase_angle).
%   |v1|  = d.
%   |v2|² = R² + d² – 2 d R cos(lat) cos(lon – phase_angle).
%   ⇒ cos(delta) = (d – R cos(lat) cos(lon – phase_angle)) / √(R² + d² – 2 d R cos(lat) cos(lon – phase_angle)).
%
% Author:  Andrea Pizzetti
% Updated: 25-Jun-2025

d_point = distanceToPointOnSphere(lon, lat, R, d, phase_angle);

cosBearing = (d - R.*cos(lat) .* cos(lon - phase_angle))./ d_point;

end