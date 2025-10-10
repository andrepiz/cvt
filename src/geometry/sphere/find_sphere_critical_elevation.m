function el = find_sphere_critical_elevation(radii, res, method)
%FIND_SPHERE_CRITICAL_ELEVATION find the worst-case elevation angle from
%the tangent plane of a set of radial locations for a dem map of given
%resolution
    
if ~exist('method','var')
    method = 'worst-case';
end

switch method
    case 'worst-case'
        a = max(radii(:));                         % maximum body radius
        b = min(radii(:));                         % minimum body radius

    case 'conservative'
        [dh_max, ix_dh_max] = max(diff(radii(:)));
        b = radii(ix_dh_max);                         % minimum body radius
        a = radii(ix_dh_max + 1);               % maximum radii difference
end

% critical elevation goes to pi/2 for resolution going to 0
ix_reverse = a.*cos(res) < b;
el = pi/2 - asin(a.*sin(res)./sqrt(a.^2 + b.^2 - 2*a.*b.*cos(res)));   
el(ix_reverse) = -el(ix_reverse);

end

