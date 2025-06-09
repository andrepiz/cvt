function [lon_subpoint, lat_subpoint] = find_direction_subpoint(vec_REF, q_REF2IAU)
% Find the longitude and latitude sub-point of a vector pointing towards
% a target in frame REF with the target oriented with respect to REF of
% given quaternion q_REF2IAU

dir_REF = vecnormalize(vec_REF);
dir_IAU = rotframe(dir_REF, q_REF2IAU);
sph = sph_coord(dir_IAU);
lon_subpoint = sph(2, :);
lat_subpoint = sph(3, :);

end