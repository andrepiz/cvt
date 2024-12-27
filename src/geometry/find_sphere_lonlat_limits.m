function [lon_lims, lat_lims, P_lims] = find_sphere_lonlat_limits(pos_c2t_TAR, dcm_TAR2CAM, R, fov, tol, flag_debug)
%FIND_SPHERE_LONLAT_LIMITS Find longitude and latitude limits covered by a
%FOV that projects onto a sphere of radius R. The FOV belongs to a camera
%placed at a certain position pos_c2t_TAR in TAR frame and oriented with
%an attitude dcm_TAR2CAM. The tolerance parameter is the difference in
%radians between following iterations under which the solution is
%considered found.

if ~exist('flag_debug','var')
    flag_debug = false;
end

if tol == 0
    error('Tolerance must be larger than 0')
end
nper = max(1, 20*ceil(1/tol));

err = tol;
lon_lims = pi/2*[-1 1];
lat_lims = pi/2*[-1 1];

while err >= tol

    % Updating with new limits
    lims_prev = [lon_lims, lat_lims];

    % Setting a finer perimeter
    nper = round(nper*1.1);

    % Computing intersections
    [P_inter_TAR, los_inter_TAR] = intersect_fov_sphere(pos_c2t_TAR, dcm_TAR2CAM, R, fov, nper);

    % Correcting number of points
    nper = size(P_inter_TAR, 2);

    % Computing limits
    sph = sph_coord(P_inter_TAR);
    [lon_lims(1), ixs_lims(1)] = min(sph(2,:));  % minimum longitude
    [lon_lims(2), ixs_lims(2)] = max(sph(2,:));  % maximum longitude
    [lat_lims(1), ixs_lims(3)] = min(sph(3,:));  % minimum latitude
    [lat_lims(2), ixs_lims(4)] = max(sph(3,:));  % maximum latitude
    
    err = max(abs([lon_lims, lat_lims] - lims_prev));
    if flag_debug
        fprintf(['\n',num2str(nper), ' samples, max error: ', num2str(rad2deg(err)), ' deg'])
    end
end

ninter = sum(~all(isnan(P_inter_TAR)));
if flag_debug
    fprintf(['\n+++ ',num2str(ninter), ' intersections found out of ', num2str(nper), '!'])
end
P_lims = P_inter_TAR(:, ixs_lims);

% Auxiliary quantities
bsDir_TAR = dcm_TAR2CAM'*[0;0;1];
d_c2t = norm(pos_c2t_TAR);
dir_c2t_TAR = pos_c2t_TAR./d_c2t;
phase_angle = atan2(-dir_c2t_TAR(2),-dir_c2t_TAR(1));
[phi_tangency, br_tangency] = find_sphere_tangent_angle(d_c2t, R);

if ninter == 0
    % No intersections found means either the FOV do not cross at all or it entirely
    % contains the sphere. To understand which situation we are, we check 
    % the angles of the FOV LOS projected on the cam-target plane with 
    % respect to boresight.

    flag_in_fov = false;
    if dot(bsDir_TAR, dir_c2t_TAR) - 1 <= eps
        % Boresight and target direction are aligned
        flag_in_fov = true;
    else
        ctPlaneDir_TAR = cross(bsDir_TAR, dir_c2t_TAR);
        losProj_TAR = project_on_plane(los_inter_TAR, ctPlaneDir_TAR);
        angBsLosProj = acos(losProj_TAR'*bsDir_TAR);  % angle between boresight and fov perimeter
        angBsT = acos(dir_c2t_TAR'*bsDir_TAR);           % angle between boresight and target    
        if angBsT <= max(angBsLosProj) + br_tangency
            flag_in_fov = true;
        end
    end
    if flag_in_fov
        lon_lims = phase_angle + phi_tangency*[-1 1];
        lat_lims = phi_tangency*[-1 1];
    end    

elseif ninter ~= nper
    % Intersect tangency circle with fov perimeter

    % Rotate to a REF frame where X is aligned with dir_t2c
    dcm_TAR2REF = euler_to_dcm([0; 0; phase_angle]);
    P_inter_REF = dcm_TAR2REF*P_inter_TAR;
    bsDir_REF = dcm_TAR2REF*bsDir_TAR;

    % Find direction of filling
    ixs_nan = all(isnan(P_inter_REF));
    %ix_first = find(~ixs_nan, 1, 'first');
    %ix_last = find(~ixs_nan, 1, 'last');
    ixs_first_nans = find(diff([0, ixs_nan])>0);
    ixs_last_nans = find(diff([ixs_nan, 0])<0);
    thetaFirst = atan2(P_inter_REF(3,ixs_first_nans(1)), P_inter_REF(2, ixs_first_nans(1)));
    thetaLast = atan2(P_inter_REF(3,ixs_last_nans(end)), P_inter_REF(2, ixs_last_nans(end)));
    thetaBs = atan2(bsDir_REF(3), bsDir_REF(2));
    if thetaLast > thetaFirst
        dir_filling = thetaBs > thetaLast | thetaBs < thetaFirst; % clockwise 
        dir_filling = +1*dir_filling -1*~dir_filling;
    else
        dir_filling = thetaBs > thetaFirst | thetaBs < thetaLast; % anticlockwise 
        dir_filling = +1*dir_filling -1*~dir_filling;
    end        

    % Fill with circles the missing points using the correct direction
    [P_inter_REF([2 3],:), ixs_nan] = fill_nans_with_arcs(P_inter_REF([2 3],:), dir_filling, flag_debug);
    P_inter_REF(1, ixs_nan) = R*cos(phi_tangency); 

    % Rotate back to TAR frame
    P_inter_TAR = dcm_TAR2REF'*P_inter_REF;

    % Computing limits
    sph = sph_coord(P_inter_TAR);
    [lon_lims(1), ixs_lims(1)] = min(sph(2,:));  % minimum longitude
    [lon_lims(2), ixs_lims(2)] = max(sph(2,:));  % maximum longitude
    [lat_lims(1), ixs_lims(3)] = min(sph(3,:));  % minimum latitude
    [lat_lims(2), ixs_lims(4)] = max(sph(3,:));  % maximum latitude
    P_lims = P_inter_TAR(:, ixs_lims);
end

if flag_debug
    k1 = R;
    k2 = R;
    k3 = R;
    nper = 100;
    [x,y,z] = sphere(100);
    figure('units','normalized','position',[0.05, 0.05, 0.8, 0.8]); 
    hold on, axis equal, camzoom(2.8), view(-pos_c2t_TAR./norm(pos_c2t_TAR))
    xlabel('X [m]'), ylabel('Y [m]'), zlabel('Z [m]')
    surf(k1*x, k1*y, k1*z, 'EdgeColor','k','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.3);
    plot3(P_inter_TAR(1,:), P_inter_TAR(2,:), P_inter_TAR(3,:),'k-','LineWidth',2)
    try
        plot3(P_inter_TAR(1, ixs_nan), P_inter_TAR(2, ixs_nan), P_inter_TAR(3, ixs_nan),'g-','LineWidth',2)
    catch
    end
    scatter3(P_lims(1,:), P_lims(2,:), P_lims(3,:), 100,'b*','LineWidth',2)
    try
        if flag_in_fov
            P_lims = cart_coord([R*ones(1,4); ...
                                 lon_lims(1), 0.5*lon_lims(1)+0.5*lon_lims(2), lon_lims(2), 0.5*lon_lims(1)+0.5*lon_lims(2); ...
                                 0.5*lat_lims(1)+0.5*lat_lims(2), lat_lims(1), 0.5*lat_lims(1)+0.5*lat_lims(2), lat_lims(2)]);
            dcm_TAR2REF = euler_to_dcm([0; 0; phase_angle]);
            Plim_plot_REF = dcm_TAR2REF*P_lims(:,1);
            P_plot_REF = zeros(3, 100);
            P_plot_REF([2, 3], :) = norm(Plim_plot_REF([2, 3],:)).*[cos(linspace(0, 2*pi, 100)); sin(linspace(0, 2*pi, 100))];
            P_plot_REF(1, :) = Plim_plot_REF(1,:);
            P_plot_TAR = dcm_TAR2REF'*P_plot_REF;
            plot3(P_plot_TAR(1, :), P_plot_TAR(2, :), P_plot_TAR(3, :),'g-','LineWidth',2)
            scatter3(P_lims(1,:), P_lims(2,:), P_lims(3,:), 100,'b*','LineWidth',2)
        end
    catch
    end
    [P_inter_TAR_plot, los_inter_TAR_plot] = intersect_fov_sphere(pos_c2t_TAR, dcm_TAR2CAM, R, fov, nper, false);
    nper = size(P_inter_TAR_plot,2);
    scatter3(P_inter_TAR_plot(1,:), P_inter_TAR_plot(2,:), P_inter_TAR_plot(3,:),'r')
    quiver3(P_inter_TAR_plot(1,:), P_inter_TAR_plot(2,:), P_inter_TAR_plot(3,:), ...
        k2*los_inter_TAR_plot(1, :), k2*los_inter_TAR_plot(2, :), k2*los_inter_TAR_plot(3, :),'r','LineWidth',1,'MarkerSize',3)
    quiver3(repmat(-pos_c2t_TAR(1,:), 1, nper), repmat(-pos_c2t_TAR(2,:), 1, nper), repmat(-pos_c2t_TAR(3,:), 1, nper), ...
       k3*los_inter_TAR_plot(1, :), k3*los_inter_TAR_plot(2, :), k3*los_inter_TAR_plot(3, :),'r','LineWidth',1,'MarkerSize',3)
end

end

