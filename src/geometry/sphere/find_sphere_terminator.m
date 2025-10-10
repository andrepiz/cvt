function [ixs, ixs_before, ixs_after] = find_sphere_terminator(lon, lat, span, k, method, tangency_angle)
%FIND_SPHERE_TERMINATOR Find the longitude and latitude sectors that belongs to
%the terminator of the sphere and within a given span angle. The parameter k tunes the reduction
%in the terminator size when approaching the poles. The tangency angle is
%used to correct the terminator regions when the light is close.

if nargin < 5
    method = 'extended';
end
if nargin < 6
    tangency_angle = pi/2;
    method = 'extended';
end
[ixs, ixs_before, ixs_after] = find_sphere_slice(lon, lat, 0, span, k, method, tangency_angle);

end