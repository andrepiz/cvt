function A = ellipsoidFrontalArea(bodyRadii, viewingDir)
% Computes projected area for (possibly repeated) ellipsoids/spheres.
% bodyRadii: (1xN, 2xN, or 3xN), or (1, 2, or 3)x1 for broadcasting
% viewingDir: 3xN (one direction vector per case)
% Returns A: 1xN areas

    N = size(viewingDir,2);
    % Expand bodyRadii to Nx needed if scalar or single case
    if size(bodyRadii,2)==1
        bodyRadii = repmat(bodyRadii,1,N);
    end
    % Expand to 3xN if needed (handle spheres and spheroids)
    if size(bodyRadii,1) == 1
        a = bodyRadii;  b = bodyRadii;  c = bodyRadii;
    elseif size(bodyRadii,1) == 2
        a = bodyRadii(1,:); b = bodyRadii(2,:); c = bodyRadii(2,:);
    else
        a = bodyRadii(1,:); b = bodyRadii(2,:); c = bodyRadii(3,:);
    end

    % Normalize viewing directions
    v = viewingDir ./ vecnorm(viewingDir);

    Q = (v(1,:).^2 ./ a.^2) + (v(2,:).^2 ./ b.^2) + (v(3,:).^2 ./ c.^2);
    A = pi .* a .* b .* c .* sqrt(Q);
end
