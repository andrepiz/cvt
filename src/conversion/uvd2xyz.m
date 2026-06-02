function [xyz, R] = uvd2xyz(K, u, v, d, sigma_u, sigma_v, sigma_d)
% Return 3D vector from pixel (u,v) + depth d.
% R = 3x3 covariance from sigmas (optional).
% K: 3x3 matrix OR 1x5 [fu fv cu cv s] (skew=0 if absent)

if numel(K) == 5
    fu = K(1); fv = K(2); cu = K(3); cv = K(4); s = K(5);
elseif numel(K) == 4
    fu = K(1); fv = K(2); cu = K(3); cv = K(4); s = 0;
else
    fu = K(1,1);  fv = K(2,2); cu = K(1,3); cv = K(2,3); s = K(1,2);
end

% Normalized coords (skew correction)
un = ((u - cu) - s*(v - cv)/fv) / fu;
vn = (v - cv) / fv;

if any(un.^2 + vn.^2 >= 1)
    error('U,V outside unit disk or K invalid');
end

zn = ones(size(un));
rn = sqrt(un.^2 + vn.^2 + zn.^2);
xyz = d .* [un; vn; zn] ./ rn;

if nargin < 5 || isempty(sigma_u) || isempty(sigma_v) || isempty(sigma_d)
    R = [];
    return
end

% Normalized sigmas (incl. skew effect on sigma_un)
sigma_un = sqrt( (sigma_u/fu)^2 + (s*sigma_v/(fu*fv))^2 );
sigma_vn = sigma_v / fv;

% Jacobian d(xyz)/d(un,vn,d)
J = [ d*(1/rn - un^2/rn^3),        -d*un*vn/rn^3,       un/rn;
     -d*un*vn/rn^3,                 d*(1/rn - vn^2/rn^3), vn/rn;
     -d*un/rn^3,                   -d*vn/rn^3,           1/rn ];

Ruvd = diag([sigma_un.^2, sigma_vn.^2, sigma_d.^2]);
R = J * Ruvd * J';
end
