function [u, v, a, b, phi] = conic2ellipse(C)
%CONIC2ELLIPSE Convert a nondegenerate image conic into ellipse parameters.
%
% C defines [u;v;1]' * C * [u;v;1] = 0.
%
% Outputs:
%   u, v - ellipse centre
%   a    - major semiaxis, a >= b
%   b    - minor semiaxis
%   phi  - major-axis orientation, modulo pi

arguments
    C (3,3) double
end

% Eliminate small asymmetries introduced numerically
C = (C + C.') / 2;

A = C(1:2,1:2);
d = C(1:2,3);

if rcond(A) < 1e-12
    error('The conic has no finite, numerically stable centre.');
end

% Ellipse centre
centre = -(A \ d);
u = centre(1);
v = centre(2);

% Conic value at the centre
centreHom = [centre; 1];
val_c = centreHom.' * C * centreHom;

if abs(val_c) <= eps(norm(C, 'fro'))
    error('The conic is degenerate or poorly conditioned.');
end

% Centred ellipse: dx' * M * dx = 1
M = -A / val_c;
M = (M + M.') / 2;

% Symmetric eigenvalue decomposition
[V, D] = eig(M, 'vector');

if any(~isfinite(D)) || any(D <= 0)
    error('The supplied conic does not represent a real ellipse.');
end

% Smallest eigenvalue corresponds to the major semiaxis
[D, idx] = sort(D, 'ascend');
V = V(:,idx);

a = 1 / sqrt(D(1));
b = 1 / sqrt(D(2));

% Major-axis orientation
majorDirection = V(:,1);
phi = atan2(majorDirection(2), majorDirection(1));

% Orientation is only defined modulo pi
phi = mod(phi + pi/2, pi) - pi/2;

% A circle has no uniquely defined orientation
if abs(a - b) <= 1e-10 * max(a,b)
    phi = 0;
end

end