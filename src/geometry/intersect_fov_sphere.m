function [P_TAR, los_TAR] = intersect_fov_sphere(pos_c2t_TAR, dcm_TAR2CAM, R, fov, nper, flag_debug)
%INTERSECT_FOV_SPHERE Find the coordinates of intersection of a FOV
%perimeter with a sphere of radius R and a camera placed at a certain
%position and attitude with respect to the target frame TAR. 
%The FOV is the full-span angle of the horizontal and vertical directions.
%If a single value is provided, the FOV is a cone, if two values are
%provided, the FOV is a frustrum.
%Sphere origin is assumed in [0; 0; 0] of target frame TAR.
%FOV Perimeter is sampled with a given number of points nper. 

if ~exist('flag_debug','var')
    flag_debug = false;
end

pos_c2t_TAR = reshape(pos_c2t_TAR, 3, 1);

% Boresight direction
bs_TAR = pos_c2t_TAR./vecnorm(pos_c2t_TAR);

% Line of sight of each point along the FOV perimeter
los_CAM = get_fov_perimeter(fov, nper, flag_debug);
los_TAR = dcm_TAR2CAM'*los_CAM;

% First point is target-to-camera position
P1_TAR = -pos_c2t_TAR;
% Second points are given by P1 plus a distance from camera along each los equal to distance of the target
P2_TAR = P1_TAR + norm(P1_TAR)*los_TAR; 
% Intersections are points in the target frame
P_TAR = intersect_line_sphere(P1_TAR, P2_TAR, [0; 0; 0], R);

% Exclude points whose camera-to-point directions have angles wrt boresight larger than
% 90 degrees
pos_c2p_TAR = pos_c2t_TAR + P_TAR;
PDir_TAR = pos_c2p_TAR./vecnorm(pos_c2p_TAR);
ixs_opposite = acos(dot(PDir_TAR, repmat(bs_TAR, 1, size(PDir_TAR, 2)))) > pi/2;
P_TAR(:, ixs_opposite) = nan;

if flag_debug
    [x, y, z] = sphere(1e2);    
    figure(); grid on, hold on, axis equal, view(P1_TAR./norm(P1_TAR)), camzoom(3)
    xlabel('x'), ylabel('y'), zlabel('z')
    surf(R*x, R*y, R*z, 'EdgeColor','k','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.2);
    scatter3(P_TAR(1,:), P_TAR(2,:), P_TAR(3,:),'r')
    quiver3(P_TAR(1,:), P_TAR(2,:), P_TAR(3,:), ...
           R*los_TAR(1, :), R*los_TAR(2, :), R*los_TAR(3, :),'r','LineWidth',1,'MarkerSize',1)
end

end

