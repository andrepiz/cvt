function [P, u] = get_sphere_tangency_circle(P1, PC, R, np)
%GET_SPHERE_TANGENTS Find the points of tangency in a sphere of radius R 
%placed at center PC with respect to a point P1. P1, PC and P are defined
%with respect to the same frame.
% P1: Point outside the sphere (3x1 vector)
% PC: Center of the sphere (3x1 vector)
% R: Radius of the sphere (scalar)
% np: Number of sampling points

pos_sph2cam_REF = PC - P1;
d = norm(pos_sph2cam_REF);
if d <= R
    P = nan(3, np);
    u = P;
    return
end

[~, ang_bearing] = find_sphere_tangent_angle(d, R);
u_CAM = get_fov_perimeter(ang_bearing*2, np, false);

% Construct CAM frame
dir_sph2cam_REF = pos_sph2cam_REF/d;
if dir_sph2cam_REF(1) ~= 0
    vecperp_REF = cross(dir_sph2cam_REF, [1; 0; 0]);
else
    vecperp_REF = cross(dir_sph2cam_REF, [0; 1; 0]);
end

zCAM_REF = dir_sph2cam_REF;
yCAM_REF = vecperp_REF./norm(vecperp_REF);
xCAM_REF = cross(yCAM_REF, zCAM_REF);
dcm_CAM2REF = [xCAM_REF, yCAM_REF, zCAM_REF];
u = dcm_CAM2REF*u_CAM;
s = sqrt(d^2 - R^2);
P = P1 + s*u;

end