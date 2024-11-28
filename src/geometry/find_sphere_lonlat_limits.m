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
lims = pi/2*[-1 1 -1 1];

while err >= tol

    lims_prev = lims;

    [P_TAR, los_TAR] = intersect_fov_sphere(pos_c2t_TAR, dcm_TAR2CAM, R, fov, nper);
    
    sph = sph_coord(P_TAR);

    [lims(1), ix_lon_min] = min(sph(2,:));
    [lims(2), ix_lon_max] = max(sph(2,:));
    [lims(3), ix_lat_min] = min(sph(3,:));
    [lims(4), ix_lat_max] = max(sph(3,:));
    
    err = max(abs(lims - lims_prev));
    if flag_debug
        disp([num2str(nper), ' points, max error: ', num2str(rad2deg(err)), ' deg'])
    end
    nper = round(nper*1.1);
end
    
lon_lims = [lims(1), lims(2)];
lat_lims = [lims(3), lims(4)];
P_lims = P_TAR(:, [ix_lon_min, ix_lon_max, ix_lat_min, ix_lat_max]);

if all(isnan([lon_lims, lat_lims]))
    % All nan means either the FOV do not cross at all or it entirely
    % contains the sphere. To understand which situation we are, we check 
    % the angles of the FOV LOS projected on the cam-target plane with 
    % respect to boresight.
    flag_in_fov = false;
    bsDir_TAR = dcm_TAR2CAM'*[0;0;1];
    d_c2t_TAR = norm(pos_c2t_TAR);
    tDir_TAR = pos_c2t_TAR./d_c2t_TAR;
    if dot(bsDir_TAR, tDir_TAR)
        % Boresight and target direction are aligned
        flag_in_fov = true;
    else
        ctPlaneDir_TAR = cross(bsDir_TAR, tDir_TAR);
        losProj_TAR = project_on_plane(los_TAR, ctPlaneDir_TAR);
        angBsLosProj = acos(losProj_TAR'*bsDir_TAR);  % angle between boresight and fov perimeter
        angBsT = acos(tDir_TAR'*bsDir_TAR);           % angle between boresight and target    sizeT = atan(R/d_c2t_TAR); % angular size of target (radius)
        sizeT = atan(R/d_c2t_TAR);                    % angular size of target radius
        if angBsT <= max(angBsLosProj) + sizeT
            flag_in_fov = true;
        end
    end
    if flag_in_fov
        angle_tangency_body = find_sphere_tangent_angle(norm(pos_c2t_TAR), R);
        lon_lims = angle_tangency_body*[-1 1];
        lat_lims = angle_tangency_body*[-1 1];
    end             
end

end
