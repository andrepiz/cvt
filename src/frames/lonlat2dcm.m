function dcm_body2point = lonlat2dcm(lon, lat)
% LONLAT2DCM(lon, lat) computes the DCM from body-fixed reference frame 
% (REF, e.g., IAU body frame) to pointing frame (e.g., CSF), where 
% x-axis points to the unit sphere surface at longitude lon (rad), 
% latitude lat (rad). Vectorized over trailing dimension of lon/lat.
%
% Inputs:
%   lon, lat - Mx1 or 1xMx... arrays (rad)
%
% Output:
%   dcm_body2point - 3x3xM... DCM array
%

x = vecnormalize([cos(lat).*cos(lon); cos(lat).*sin(lon); sin(lat)]);
z_init = vecnormalize([sin(lon); zeros(size(lon)); cos(lon)]);  
y = cross(z_init, x);
ixs_sing = vecnorm(y) <= eps;
if any(ixs_sing)
    y(:, ixs_sing) = cross(repmat([-1; 0; 0], 1, size(z_init, 2)), x);
end
y = vecnormalize(y);
z = vecnormalize(cross(x, y));
dcm_body2point = zeros(3,3,size(x,2));

dcm_body2point(1,1:3,:) = x;
dcm_body2point(2,1:3,:) = y;
dcm_body2point(3,1:3,:) = z;

end
