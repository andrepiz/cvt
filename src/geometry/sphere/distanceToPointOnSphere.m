function d_point = distanceToPointOnSphere(lon, lat, R, d, phase_angle)
% DISTANCETOPOINTONSPHERE  Find the distance between a point at longitude 
% and latitude coordinates on a sphere and an observer placed at a given 
% phase angle. 
%
%   d_point = distanceToPointOnSphere(lon, lat, R, d, phase_angle) 
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
%   • The sphere is centered at the origin O = (0, 0, 0).
%   • The external point P lies in the XY-plane:
%         P = (d * cos(phase_angle), d * sin(phase_angle), 0)
%   • The surface point S on the sphere is defined in spherical coordinates:
%         S = (R * cos(lat) * cos(lon), 
%              R * cos(lat) * sin(lon), 
%              R * sin(lat))
%
% Derivation
% ----------
%   Let v1 = vector from S to P: v1 = P - S
%   Then, using the Euclidean distance and knowing cos(alpha - beta) = cos alpha cos beta + sin alpha sin beta
%
%     |P - S|^2 = |P|^2 + |S|^2 - 2 * P · S
%              = d^2 + R^2 - 2*d*R*cos(lat)*cos(lon - phase_angle)
%
%   Therefore,
%     d_point = sqrt(R^2 + d^2 - 2*d*R*cos(lat)*cos(lon - phase_angle))
%
% Author:  Andrea Pizzetti
% Updated: 25-Jun-2025
    
% proof: @lon = phase_angle, lat = tangency_angle = pi/2 - asin(R/d)
% --> d_point = sqrt(R^2 + d^2 - 2*R*d*cos(pi/2-asin(R/d)) = sqrt(d^2 - R^2) = d_inter

d_point = sqrt( R.^2 + d.^2 - 2*R.*d.*cos(lat).*cos(lon - phase_angle) );

end
