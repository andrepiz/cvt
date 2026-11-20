function conicMat = ellipsoid2conicMat_IAU(pos_cam2body_IAU, dcm_CAM2IAU, K, R)
%ELLIPSOID2CONICMAT_IAU Conic of an ellipsoid, REF = IAU (body-principal).
%
% conicMat = ellipsoid2conicMat_IAU(pos_cam2body_IAU, dcm_CAM2IAU, K, R)
%
% Inputs
%   pos_cam2body_IAU - 3x1 camera-to-body-center vector, in IAU frame.
%   dcm_CAM2IAU      - 3x3 DCM from CAM to IAU.
%   K                - 3x3 camera intrinsics.
%   R                - 1x1/1x2/1x3 semiaxes.
%
% Output
%   conicMat         - 3x3 conic in pixel coordinates (primal conic).

% Normalize R
nR = numel(R);
if nR == 1
    aP = R; bP = R; cP = R;
elseif nR == 2
    aP = R(1); bP = R(1); cP = R(2);
elseif nR == 3
    aP = R(1); bP = R(2); cP = R(3);
else
    error('R must have 1, 2, or 3 elements.');
end

% Shape matrix in IAU (principal) frame: x' Ap x = 1
Ap = diag([1/aP^2, 1/bP^2, 1/cP^2]);

% DCMs
C_CI = dcm_CAM2IAU;        % CAM -> IAU
C_IC = C_CI.';             % IAU -> CAM

% Camera-to-body vector in IAU and in CAM
r_I = pos_cam2body_IAU;    % 3x1
r_C = C_IC * r_I;          % 3x1

% Shape matrix in camera frame: Ac = C_IC' * Ap * C_IC
Ac = C_IC.' * Ap * C_IC;   % = C_CI * Ap * C_CI'

% Ensure we are outside the ellipsoid: r' Ac r > 1
check = r_C.' * Ac * r_C;
if check <= 1
    warning('Camera inside or on the ellipsoid (r''Ac r <= 1). Conic may be degenerate.');
end

% Dual conic in normalized image coordinates (up to scale):
% C*_norm ∝ Ac - (Ac r)(Ac r)' / (r' Ac r - 1)
den = check - 1;
if abs(den) < eps
    warning('Degenerate configuration: r''Ac r ≈ 1.');
    den = eps;
end
Ac_r = Ac * r_C;
Cstar_norm = Ac - (Ac_r * Ac_r.') / den;  % 3x3, symmetric

% Convert dual conic in normalized coords to primal conic in pixels:
% C ∝ K^{-T} * inv(C*_norm) * K^{-1}
Cstar_norm = (Cstar_norm + Cstar_norm.') / 2; % enforce symmetry
Cstar_inv = inv(Cstar_norm);                  % primal in normalized coords

Kinv = inv(K);
conicMat = Kinv.' * Cstar_inv * Kinv;
conicMat = (conicMat + conicMat.') / 2;

% Normalize
if abs(conicMat(3,3)) > eps
    conicMat = conicMat / conicMat(3,3);
end

end