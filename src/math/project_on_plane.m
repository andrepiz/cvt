function vec_proj = project_on_plane(vec, dir)
% Project vector vec onto a plane with direction dir

dir = reshape(dir, 3, []);
nd = size(dir, 2);

vec = reshape(vec, 3, []);
nv = size(vec, 2);

dir = dir./vecnorm(dir);

if nd == nv
    vec_proj = vec - (dot(dir, vec)).*dir;
else
    error('Dimensions of vec and dir are not consistent')
end

end