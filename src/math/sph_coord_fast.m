function sph_coord = sph_coord_fast(x)
% Convert Cartesian coordinates to spherical coordinates
% Input: x (3 x N) array
if size(x,1)~=3
    error('Input must have 3 rows')
end

sph_coord = zeros(size(x));

sph_coord(2,:) = atan2(x(2,:), x(1,:));
pos_XY = hypot(x(1,:), x(2,:));
sph_coord(1,:) = hypot(pos_XY, x(3,:));
sph_coord(3,:) = atan2(x(3,:), pos_XY);

end
