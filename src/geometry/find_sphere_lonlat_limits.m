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

nper = 20*round(1/tol);

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
phase_angle = acos(-dir_c2t_TAR(1));

if ninter == 0
    % No intersections found means either the FOV do not cross at all or it entirely
    % contains the sphere. To understand which situation we are, we check 
    % the angles of the FOV LOS projected on the cam-target plane with 
    % respect to boresight.

    flag_in_fov = false;
    if dot(bsDir_TAR, dir_c2t_TAR)
        % Boresight and target direction are aligned
        flag_in_fov = true;
    else
        ctPlaneDir_TAR = cross(bsDir_TAR, dir_c2t_TAR);
        losProj_TAR = project_on_plane(los_inter_TAR, ctPlaneDir_TAR);
        angBsLosProj = acos(losProj_TAR'*bsDir_TAR);  % angle between boresight and fov perimeter
        angBsT = acos(dir_c2t_TAR'*bsDir_TAR);           % angle between boresight and target    sizeT = atan(R/d_c2t_TAR); % angular size of target (radius)
        sizeT = atan(R/d_c2t);                    % angular size of target radius
        if angBsT <= max(angBsLosProj) + sizeT
            flag_in_fov = true;
        end
    end
    if flag_in_fov
        angle_tangency_body = find_sphere_tangent_angle(d_c2t, R);
        lon_lims = phase_angle + angle_tangency_body*[-1 1];
        lat_lims = angle_tangency_body*[-1 1];
    end    

elseif ninter ~= nper
    % Intersect tangency circle with fov perimeter

    % Rotate to a REF frame where X is aligned with dir_t2c
    dcm_TAR2REF = euler_to_dcm([0; 0; phase_angle]);
    P_inter_REF = dcm_TAR2REF*P_inter_TAR;

    % Fill with circles the missing points
    [P_inter_REF([2 3],:), ixs_filled] = fill_nans_with_arcs(P_inter_REF([2 3],:), 1, flag_debug);
    angle_tangency_body = find_sphere_tangent_angle(norm(pos_c2t_TAR), R);
    P_inter_REF(1, ixs_filled) = R*cos(angle_tangency_body); 

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
    k2 = R/2;
    k3 = R/2;
    nper = 100;
    [x,y,z] = sphere(100);
    figure(); grid on, hold on, axis equal, camzoom(1), view(-bsDir_TAR)
    surf(k1*x, k1*y, k1*z, 'EdgeColor','k','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.2);
    plot3(P_inter_TAR(1,:), P_inter_TAR(2,:), P_inter_TAR(3,:),'k-')
    scatter3(P_lims(1,:), P_lims(2,:), P_lims(3,:),'b*','LineWidth',5)
    [P_inter_TAR_plot, los_inter_TAR_plot] = intersect_fov_sphere(pos_c2t_TAR, dcm_TAR2CAM, R, fov, nper, false);
    scatter3(P_inter_TAR_plot(1,:), P_inter_TAR_plot(2,:), P_inter_TAR_plot(3,:),'r')
    quiver3(P_inter_TAR_plot(1,:), P_inter_TAR_plot(2,:), P_inter_TAR_plot(3,:), ...
        k2*los_inter_TAR_plot(1, :), k2*los_inter_TAR_plot(2, :), k2*los_inter_TAR_plot(3, :),'r','LineWidth',1,'MarkerSize',1)
    quiver3(repmat(-pos_c2t_TAR(1,:), 1, nper), repmat(-pos_c2t_TAR(2,:), 1, nper), repmat(-pos_c2t_TAR(3,:), 1, nper), ...
       k3*los_inter_TAR_plot(1, :), k3*los_inter_TAR_plot(2, :), k3*los_inter_TAR_plot(3, :),'r','LineWidth',1,'MarkerSize',1)
end

end

