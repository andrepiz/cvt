function [ixs, ixs_before, ixs_after] = find_sphere_slice(lon, lat, phase_angle, span, k, method)
%FIND_SPHERE_SLICE Find the longitude and latitude sectors that belongs to
%the slice of a sphere at a given phase angle with respect to the zero longitude
%point. The slice points are extracted within a given span angle with respect
%to the slice meridinan. The parameter k tunes the reduction
%in the slice size when approaching the poles.

if nargin < 6
    method = 'extended';
end

% Wrap longitude points to the slice meridian
lon = wrapToPi(lon-phase_angle);

% Extract only points around slice
switch method
    case 'simple'
        % Constrain the angle between the slice line to be under a
        % certain threshold given by span angle
        thr = abs(span).*cos(lat).^k;
        ixs = abs(abs(lon) - pi/2) < thr;
        ixs_before = ~ixs & abs(lon) < pi/2;
        ixs_after = ~ixs_before & ~ixs; 
    case 'extended'
        % Constrain the longitude value to fall between maximum span at
        % equator and any value at the poles
        thr = (pi/2 - abs(span)).*cos(lat).^k;
        ixs = min(abs(lon), abs(abs(lon) - pi)) > thr;
        ixs_before = ~ixs & abs(lon) < pi/2;
        ixs_after = ~ixs_before & ~ixs; 
    otherwise
        error('Method not supported')

end

end