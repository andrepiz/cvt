function [ixs, ixs_before, ixs_after] = find_sphere_terminator(lon, lat, span, k, method)
%FIND_SPHERE_TERMINATOR Find the longitude and latitude sectors that belongs to
%the terminator of the sphere and within a given span angle. The parameter k tunes the reduction
%in the terminator size when approaching the poles.

[ixs, ixs_before, ixs_after] = find_sphere_slice(lon, lat, 0, span, k, method);

end