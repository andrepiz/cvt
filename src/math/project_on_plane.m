function vec_proj = project_on_plane(vec, dir)
% Project vector vec onto a plane with direction dir

dir = reshape(dir, 3, []);
nd = size(dir, 2);

vec = reshape(vec, 3, []);
nv = size(vec, 2);

dir = dir./vecnorm(dir);

if nd == nv
    vec_proj = vec - (dot(dir, vec)).*dir;
elseif nd == 1
    vec_proj = vec - dot(repmat(dir, 1, nv), vec).*repmat(dir, 1, nv);
elseif nv == 1
    vec_proj = repmat(vec, 1, nv) - dot(dir, repmat(vec, 1, nv)).*dir;
else
    error('Dimensions of vec and dir are not consistent')
end

end