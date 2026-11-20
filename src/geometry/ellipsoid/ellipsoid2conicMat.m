function conicMat = ellipsoid2conicMat(pos_origin2body_REF, pos_body2cam_REF, dcm_REF2CAM, K, R)
%ELLIPSOID2CONICMAT Compute image conic of an ellipsoid (sphere/spheroid/triaxial).
%
% conicMat = ellipsoid2conicMat(pos_origin2body_REF, pos_body2cam_REF, dcm_REF2CAM, K, R)
%
% Inputs
%   pos_origin2body_REF - 3x1 vector from reference frame origin to body center,
%                         expressed in the reference frame [xc; yc; zc].
%   pos_body2cam_REF    - 3x1 vector from body center to camera, expressed in REF.
%   dcm_REF2CAM         - 3x3 DCM from reference frame to camera frame.
%   K                   - 3x3 camera intrinsic matrix.
%   R                   - 1x1 (sphere), 1x2 (spheroid), or 1x3 (triaxial) semiaxes.
%                         1 elem: a = b = c = R
%                         2 elem: a = b = R(1), c = R(2)
%                         3 elem: a = R(1), b = R(2), c = R(3)
%
% Output
%   conicMat            - 3x3 symmetric conic matrix in pixel coordinates,
%                         such that p' * conicMat * p = 0 for p = [u; v; 1].

% Body center in REF
rO = pos_origin2body_REF;          % 3x1
pos_origin2cam_REF = pos_origin2body_REF + pos_body2cam_REF;

% Normalize R to [a, b, c]
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

% Ellipsoid shape matrix in body principal frame:
% x_p' * Ap * x_p = 1, with Ap = diag(1/a^2, 1/b^2, 1/c^2)
Ap = diag([1/aP^2, 1/bP^2, 1/cP^2]);

% Assume principal axes aligned with REF; if not, insert rotation here:
% Cbp = dcm_REF2BODY;  % REF -> body principal
% A_REF = Cbp.' * Ap * Cbp;
A_REF = Ap;

% Quadric in REF homogeneous coordinates (4x4):
% Q_REF = [ A_REF,  -A_REF * rO;
%          -rO' * A_REF,  rO' * A_REF * rO - 1 ]
Aro = A_REF * rO;
roAro = rO' * Aro;                 % scalar
Q_REF = [ A_REF,       -Aro;
         -Aro.',  roAro - 1 ];     % 4x4 symmetric

% Camera pose: from REF to CAM
tSC = -dcm_REF2CAM * pos_body2cam_REF;  % camera center in CAM frame
rotSC = dcm_REF2CAM;                    % REF -> CAM

% Camera projection matrix (3x4): P = K * [R | t]
P = K * [rotSC, tSC];                   % 3x4

Cstar = P * (Q_REF \ P.');
conicMat = Cstar \ eye(3);
conicMat = (conicMat + conicMat.') / 2;

% Optional normalization: set C(3,3) = 1 if non-zero
if abs(conicMat(3,3)) > eps
    conicMat = conicMat / conicMat(3,3);
end

end