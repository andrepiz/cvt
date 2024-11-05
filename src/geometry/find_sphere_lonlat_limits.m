function [lon_lims, lat_lims, P_lims] = find_sphere_lonlat_limits(pos_c2t_TAR, dcm_TAR2CAM, R, fov, tol, flag_debug)
%FIND_SPHERE_LONLAT_LIMITS Find longitude and latitude limits covered by a
%FOV that projects onto a sphere of radius R. The FOV belongs to a camera
%placed at a certain position pos_c2t_TAR in TAR frame and oriented with
%an attitude dcm_TAR2CAM. The tolerance parameter is the difference in
%radians between following iterations under which the solution is
%considered found.
%Boresight of camera is assumed aligned to +Z direction
% bs = [0; 0; 1];

% Number of points to sample depends on how much tilted is the camera with
% respect to the direction aligned with the line linking the center of the
% sphere with the camera
% tiltAngle = norm(fov) + acos(dot(pos_c2t_TAR./norm(pos_c2t_TAR), dcm_TAR2CAM'*bs));
% % Maximum tilt angle
% maxAngle = pi/2;
% % Maximum sampling number of points 
% nperMax = 1e4;
% 
% nper = round(tiltAngle/maxAngle*nperMax);

nper = 20*round(1/tol);

err = tol;
lims = pi/2*[-1 1 -1 1];

while err >= tol

    lims_prev = lims;

    P = intersect_fov_sphere(pos_c2t_TAR, dcm_TAR2CAM, R, fov, nper);
    
    sph = sph_coord(P);

    [lims(1), ix_lon_min] = min(sph(2,:));
    [lims(2), ix_lon_max] = max(sph(2,:));
    [lims(3), ix_lat_min] = min(sph(3,:));
    [lims(4), ix_lat_max] = max(sph(3,:));
    
    err = max(abs(lims - lims_prev));

    disp([num2str(nper), ' points, max error: ', num2str(rad2deg(err)), ' deg'])

    nper = round(nper*1.1);
end
    
lon_lims = [lims(1), lims(2)];
lat_lims = [lims(3), lims(4)];
P_lims = P(:, [ix_lon_min, ix_lon_max, ix_lat_min, ix_lat_max]);

end
