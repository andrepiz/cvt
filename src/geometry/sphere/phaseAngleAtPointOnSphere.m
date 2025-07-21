function phase_angle_point = phaseAngleAtPointOnSphere(lon, lat, R, d, phase_angle)
% PHASEANGLEATPOINTONSPHERE  Local phase angle on a spherical surface (radians).
%
%   phase_angle_local = phaseAngleAtPointOnSphere(lon, lat, R, d, phase_angle)
%
%   Inputs
%   ------
%   lon, lat   : planetocentric longitude & latitude of the surface point
%                (radians; 0 lon = +X, +lat = +Z).
%   R          : radius of the sphere (same units as d).
%   d          : observer distance from the sphere’s centre.
%   phase_angle: planetary phase angle α = ∠(Sun‑centre‑observer) [radians].
%                The observer lies in the equatorial (XY) plane at this
%                azimuth measured from +X (sub‑solar) toward +Y.
%
%   Output
%   -------
%   phase_angle_local : local phase angle ϕ at (lon, lat) in radians.
%
%   Geometry
%   --------
%       Sun unit vector      s = [1 0 0]
%       Observer position    o = d · [cos α, sin α, 0]
%
%   Example
%   -------
%       R_earth = 6.371e6;         % m
%       d_obs   = 4.0e8;           % 400 000 km from centre
%       alpha   = deg2rad(90);     % half‑phase
%       lon = 0;  lat = 0;         % sub‑solar point
%       phi = phaseAngleAtPointOnSphere(lon, lat, R_earth, d_obs, alpha)


d_point = distanceToPointOnSphere(lon, lat, R, d, phase_angle);

phase_angle_point = acos( (d.*cos(phase_angle) - R.*cos(lat).*cos(lon)) ./ d_point);
end
