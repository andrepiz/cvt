cvt_install()

%%

% flags
flag_scatter_pose = true;
flag_create_gif = false;
ns = 20;
tol = deg2rad(0.05);

% cam
muPixel = 5.5e-6;
f = 10e-3;
fov = deg2rad([15, 10]);
%fov = deg2rad([70, 170]);
%fov = deg2rad([50, 70]);
%fov = deg2rad([4, 10]);
%fov = deg2rad([40]);
res_px = 2/muPixel*f*tan(fov/2);

% geometry
phase_angle = 0.2;
Rbody = 1737.4;

if flag_scatter_pose
    eul_CAMI2CAM = (2*rand(3, ns)-1).*[deg2rad(20);deg2rad(20);pi];
    pos_body2cam_TAR = 2*(rand(1, ns)+1)*Rbody.*[cos(phase_angle); sin(phase_angle); 0];
else
    %eul_CAMI2CAM = deg2rad([22; 12; 20]);
    eul_CAMI2CAM = deg2rad([16; -1; 20]);
    %eul_CAMI2CAM = deg2rad([-50; -20; -3]);
    eul_CAMI2CAM = deg2rad([-10; -5; 0]);
    %eul_CAMI2CAM = deg2rad([1; -16; 50]);
    %eul_CAMI2CAM = deg2rad([22; 52; 20]);
    %eul_CAMI2CAM = deg2rad([-22; 52; -150]);
    %eul_CAMI2CAM = deg2rad([42; 12; 20]);
    pos_body2cam_TAR = 4*Rbody*[cos(phase_angle); sin(phase_angle); 0];
    %pos_body2cam_TAR = 1.1*Rbody*[cos(phase_angle); sin(phase_angle); 0];
    % 
    % %---TRICKY
    % fov = deg2rad([140, 170]);
    % eul_CAMI2CAM = deg2rad([22; 12; 20]);
    % pos_body2cam_TAR = 2*Rbody*[0.8; 0.7; 0];
    % %---
end

if flag_scatter_pose

    nsph = 100;
    [x, y, z] = sphere(nsph);

    k1 = Rbody;
    k2 = Rbody;
    k3 = Rbody/2;

    % Plot the sphere
    figure('units','normalized','position',[0.1, 0.1, 0.8, 0.9]); 
    hold on, axis equal, camzoom(5)
    xlabel('x'), ylabel('y'), zlabel('z')
    surf(k1*x, k1*y, k1*z, 'EdgeColor','k','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.3);
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    view(pos_body2cam_TAR(:,1)/norm(pos_body2cam_TAR(:,1)))
    xlim(4*Rbody*[-1, 1]), ylim(4*Rbody*[-1, 1]), zlim(4*Rbody*[-1, 1]);

    if flag_create_gif
        gif('fov_intersection_scatter_pose_2.gif','DelayTime',0.5)
    end
    for ix = 1:size(pos_body2cam_TAR, 2)
        pos_cam2body_TAR = -pos_body2cam_TAR(:, ix);
        dir_cam2body_TAR = pos_cam2body_TAR./vecnorm(pos_cam2body_TAR);
        
        zCAM_TAR = dir_cam2body_TAR;
        yCAM_TAR = -[0; 0; 1];
        xCAM_TAR = cross(yCAM_TAR, zCAM_TAR);
        dcm_TAR2CAMI = [xCAM_TAR, yCAM_TAR, zCAM_TAR]';
        dcm_CAMI2CAM = euler_to_dcm(eul_CAMI2CAM(:, ix));
        dcm_TAR2CAM = dcm_CAMI2CAM*dcm_TAR2CAMI;
        
        nper = 50;
        nsph = 100;
        
        bsP1_TAR = pos_body2cam_TAR(:, ix);        % P1 at camera
        bsDir_TAR = dcm_TAR2CAM'*[0;0;1];
        bsP2_TAR = bsP1_TAR + Rbody*bsDir_TAR; % P2 at 1 radius distance from camera along boresight
        bs_intersection_TAR = intersect_line_sphere(bsP1_TAR, bsP2_TAR, [0;0;0], Rbody);
        [fov_intersection_TAR, los_TAR] = intersect_fov_sphere(pos_cam2body_TAR, dcm_TAR2CAM, Rbody, fov, nper, false);    
        [lon_lims_TAR, lat_lims_TAR, fov_lims_TAR] = find_sphere_lonlat_limits(pos_cam2body_TAR, dcm_TAR2CAM, Rbody, fov, tol, false);
        nper = size(fov_intersection_TAR, 2);

        % PLOT
        objs{1} = scatter3(bs_intersection_TAR(1,:), bs_intersection_TAR(2,:), bs_intersection_TAR(3,:),'g','LineWidth',3);
        objs{2} = quiver3(bsP1_TAR(1,:), bsP1_TAR(2,:), bsP1_TAR(3,:), ...
           k2*bsDir_TAR(1, :), k2*bsDir_TAR(2, :), k2*bsDir_TAR(3, :),'g','LineWidth',2,'MarkerSize',3);
        objs{3} =scatter3(fov_intersection_TAR(1,:), fov_intersection_TAR(2,:), fov_intersection_TAR(3,:),'r');
        objs{4} =quiver3(fov_intersection_TAR(1,:), fov_intersection_TAR(2,:), fov_intersection_TAR(3,:), ...
           k2*los_TAR(1, :), k2*los_TAR(2, :), k2*los_TAR(3, :),'r','LineWidth',1,'MarkerSize',1);
        objs{5} =quiver3(repmat(pos_body2cam_TAR(1,ix), 1, nper), repmat(pos_body2cam_TAR(2,ix), 1, nper), repmat(pos_body2cam_TAR(3,ix), 1, nper), ...
           k3*los_TAR(1, :), k3*los_TAR(2, :), k3*los_TAR(3, :),'r','LineWidth',1,'MarkerSize',1);
        objs{6} =scatter3(fov_lims_TAR(1,:), fov_lims_TAR(2,:), fov_lims_TAR(3,:),'b*','LineWidth',5);
        
        pause(0.2)
        if flag_create_gif
            gif
        end
        if ix < size(pos_body2cam_TAR, 2)
            cellfun(@delete, objs)
        end
    
    end

else
    pos_cam2body_TAR = -pos_body2cam_TAR;
    dir_cam2body_TAR = pos_cam2body_TAR./vecnorm(pos_cam2body_TAR);
    
    zCAM_TAR = dir_cam2body_TAR;
    yCAM_TAR = -[0; 0; 1];
    xCAM_TAR = cross(yCAM_TAR, zCAM_TAR);
    dcm_TAR2CAMI = [xCAM_TAR, yCAM_TAR, zCAM_TAR]';
    dcm_CAMI2CAM = euler_to_dcm(eul_CAMI2CAM);
    dcm_TAR2CAM = dcm_CAMI2CAM*dcm_TAR2CAMI;
    
    nper = 100;
    nsph = 100;
    
    bsP1_TAR = pos_body2cam_TAR;        % P1 at camera
    bsDir_TAR = dcm_TAR2CAM'*[0;0;1];
    bsP2_TAR = bsP1_TAR + Rbody*bsDir_TAR; % P2 at 1 radius distance from camera along boresight
    bs_intersection_TAR = intersect_line_sphere(bsP1_TAR, bsP2_TAR, [0;0;0], Rbody);
    
    [fov_intersection_TAR, los_TAR] = intersect_fov_sphere(pos_cam2body_TAR, dcm_TAR2CAM, Rbody, fov, nper, true);
    
    tic
    [lon_lims_TAR, lat_lims_TAR, fov_lims_TAR] = find_sphere_lonlat_limits(pos_cam2body_TAR, dcm_TAR2CAM, Rbody, fov, tol, true);
    toc
end

