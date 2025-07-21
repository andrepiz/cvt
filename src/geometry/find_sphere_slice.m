function [ixs, ixs_before, ixs_after] = find_sphere_slice(lon, lat, phase_angle, span, k, method, tangency_angle)
%FIND_SPHERE_SLICE Find the longitude and latitude sectors that belongs to
%the slice of a sphere at a given phase angle with respect to the zero longitude
%point. The slice points are extracted within a given span angle with respect
%to the slice meridinan. The parameter k tunes the reduction
%in the slice size when approaching the poles. The tangency angle is
%used to correct the slice regions when the light is close.

if nargin < 6
    method = 'extended';
end
if nargin < 7
    tangency_angle = pi/2;
    method = 'extended';
end

% Wrap longitude points to the slice meridian
lon = wrapToPi(lon-phase_angle);

% Select only those points whose longitude difference with respect 
% to the longitude of the slice is under a certain threshold. 
switch method
    case 'simple'
        % the threshold is zero at the polar regions and equal to the span
        % angle at the equator.
        thr = abs(span).*cos(lat).^k;
        ixs = abs(abs(lon) - tangency_angle) < thr;
        ixs_before = ~ixs & abs(lon) < tangency_angle;
        ixs_after = ~ixs_before & ~ixs; 
    case 'extended'
        % the threshold is any value at the polar regions and equal to the
        % span angle at the equator
        thr = pi - (pi - abs(span)).*cos(lat).^k;
        ixs = abs(abs(lon) - tangency_angle) < thr;
        ixs_before = ~ixs & abs(lon) < tangency_angle;
        ixs_after = ~ixs_before & ~ixs; 
        % thr = (pi/2 - abs(span)).*cos(lat).^k;
        % ixs = min(abs(lon), abs(abs(lon) - pi)) > thr;
        % ixs_before = ~ixs & abs(lon) < pi/2;
        % ixs_after = ~ixs_before & ~ixs; 
    otherwise
        error('Method not supported')

end

end