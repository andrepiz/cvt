function A = ellipsoidFrontalArea(bodyRadii, viewingDir)
% Computes projected area for (possibly repeated) ellipsoids/spheres.
% bodyRadii: 3x1 or 1x3
% viewingDir: 3xN (one direction vector per case) vector from center of the
% body to the observer
% Returns A: 1xN areas

N = size(viewingDir,2);

bodyRadii = reshape(bodyRadii, 1, []);
if numel(bodyRadii) == 1
    a = bodyRadii;  b = bodyRadii;  c = bodyRadii;
elseif numel(bodyRadii) == 2
    a = bodyRadii(1); b = bodyRadii(1); c = bodyRadii(2);
elseif numel(bodyRadii) == 3
    a = bodyRadii(1); b = bodyRadii(2); c = bodyRadii(3);
else
    error('bodyRadii should be a vector of 1, 2 or 3 elements')
end

% Normalize viewing directions
v = viewingDir ./ vecnorm(viewingDir);

Q = (v(1,:).^2 ./ a.^2) + (v(2,:).^2 ./ b.^2) + (v(3,:).^2 ./ c.^2);
A = pi .* a .* b .* c .* sqrt(Q);

end
