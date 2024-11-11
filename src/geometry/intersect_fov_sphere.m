function [P, los] = intersect_fov_sphere(pos_c2t_TAR, dcm_TAR2CAM, R, fov, nper, flag_debug)
%INTERSECT_FOV_SPHERE Find the coordinates of intersection of a FOV
%perimeter with a sphere of radius R and a camera placed at a certain
%position and attitude with respect to the target frame TAR. 
%Sphere origin is assumed in [0; 0; 0] of target frame TAR.
%FOV Perimeter is sampled with a given number of points nper. 

if ~exist('flag_debug','var')
    flag_debug = false;
end

pos_c2t_TAR = reshape(pos_c2t_TAR, 3, 1);

los_CAM = get_fov_perimeter(fov, nper, flag_debug);

los = dcm_TAR2CAM'*los_CAM;

losP1 = -pos_c2t_TAR;
losP2 = losP1 + norm(losP1)*los; % P2 at distance from camera along boresight equal to distance of the target 

P = intersect_line_sphere(losP1, losP2, [0; 0; 0], R);

if flag_debug
    [x, y, z] = sphere(1e2);    
    figure(); grid on, hold on, axis equal, view(losP1./norm(losP1)), camzoom(3)
    xlabel('x'), ylabel('y'), zlabel('z')
    surf(R*x, R*y, R*z, 'EdgeColor','k','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.2);
    scatter3(P(1,:), P(2,:), P(3,:),'r')
    quiver3(P(1,:), P(2,:), P(3,:), ...
           R*los(1, :), R*los(2, :), R*los(3, :),'r','LineWidth',1,'MarkerSize',1)
end

end

