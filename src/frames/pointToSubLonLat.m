function [q_IAU2CAM, dcm_IAU2CAM] = pointToSubLonLat(subLon, subLat)
%POINTTOSUBLONLAT Set camera so to point a sub-longitude and sub-latitude
%point at nadir. The camera vertical coordinates axes always point the south pole.
%When sub-latitude is 90 deg or -90 deg, singularity is solved aligning 
%the camera axes to IAU.

numTol = 1e-11;
n = size(subLon, 2);

z = -vecnormalize([cos(subLat).*cos(subLon); cos(subLat).*sin(subLon); sin(subLat)]);
y_init = repmat([0; 0; -1], 1, n);  
x = cross(y_init, z);
ixs_sing = abs(dot(y_init, z) - 1) <= numTol;
if any(ixs_sing)
    x(:, ixs_sing) = repmat([1; 0; 0], 1, n);
end
x = vecnormalize(x);
y = vecnormalize(cross(z, x));
dcm_IAU2CAM = zeros(3, 3, n);

dcm_IAU2CAM(1, 1:3, :) = x;
dcm_IAU2CAM(2, 1:3, :) = y;
dcm_IAU2CAM(3, 1:3, :) = z;

q_IAU2CAM = dcm_to_quat(dcm_IAU2CAM);

end

