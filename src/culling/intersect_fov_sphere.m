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

% adimensionalize for numerical stability
adim = R;

posAdim_c2t_TAR = reshape(pos_c2t_TAR/adim, 3, 1);
dAdim_c2t_TAR = norm(posAdim_c2t_TAR);
if dAdim_c2t_TAR - 1 < -eps
     error('Camera is inside the sphere')
end

% Boresight direction
bs_TAR = dcm_TAR2CAM'*[0; 0; 1];

% Line of sight of each point along the FOV perimeter
los_CAM = get_fov_perimeter(fov, nper, flag_debug);
los_TAR = dcm_TAR2CAM'*los_CAM;

% First point is target-to-camera position
PAdim1_TAR = -posAdim_c2t_TAR;
% Second points are given by P1 plus a distance from camera along each los equal to distance of the target
PAdim2_TAR = PAdim1_TAR + dAdim_c2t_TAR*los_TAR; 
% Intersections are points in the target frame
PAdim_TAR = intersect_line_sphere(PAdim1_TAR, PAdim2_TAR, [0; 0; 0], 1);

% Exclude points whose camera-to-point directions have angles wrt boresight larger than
% 90 degrees
pos_c2p_TAR = posAdim_c2t_TAR + PAdim_TAR;
PDir_TAR = pos_c2p_TAR./vecnorm(pos_c2p_TAR);
ixs_opposite = bs_TAR'*PDir_TAR < 0;
PAdim_TAR(:, ixs_opposite) = nan;

% Dimensionalize
P_TAR = PAdim_TAR*adim;

if flag_debug
    [x, y, z] = sphere(1e2);    
    figure(); grid on, hold on, axis equal, view(PAdim1_TAR./norm(PAdim1_TAR)), camzoom(3)
    xlabel('x'), ylabel('y'), zlabel('z')
    surf(x, y, z, 'EdgeColor','k','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.2);
    scatter3(PAdim2_TAR(1,:), PAdim2_TAR(2,:), PAdim2_TAR(3,:),'b')
    scatter3(PAdim_TAR(1,:), PAdim_TAR(2,:), PAdim_TAR(3,:),'r')
    quiver3(PAdim_TAR(1,:), PAdim_TAR(2,:), PAdim_TAR(3,:), ...
           los_TAR(1, :), los_TAR(2, :), los_TAR(3, :),'r','LineWidth',1,'MarkerSize',1)
    cameratoolbar('show')
end

end

