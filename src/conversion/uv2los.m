function los = uv2los(K, u, v)
% Return the 3D los vector from the camera pixel coordinates U, V

los = uvd2xyz(K, u, v, ones(1, length(u)));

end