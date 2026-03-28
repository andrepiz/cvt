%%
cvt_install()

%% Inputs
AU = 149597870707;
Rbody = 1737.4e3;
phase_angle = deg2rad(5);
np = 1e6;   % sectors
nrs = 10;  % ray sampling points
pos_body2star_CSF = AU*[1;0;0];
pos_body2cam_CSF = 80e6*[cos(phase_angle);sin(phase_angle);0];
dcm_CSF2IAU = euler_to_dcm([pi/4; pi/4; pi/4]);
%dcm_CSF2IAU = euler_to_dcm([0; 0; 0]);
filename = 'ldem_4.tif';
method_intersection = 'sampling';
method_sampling = 'linspace';
nworkers = 12;

method_threshold = 'extended';
k_threshold = 0.2;

% 3d plots
plot_mask_only = true;
cmap = 'parula';
sm = 2;

%% Prepare data

np = round(sqrt(np)).^2;
dir_body2cam_CSF = pos_body2cam_CSF./vecnorm(pos_body2cam_CSF);
dir_body2star_CSF = pos_body2star_CSF./vecnorm(pos_body2star_CSF);

% Using Moon radius and DEM
img = imread(filename);
map_img = digital2analog(img, 1, 1, 1);
map_img = matop(map_img, 'scale', 1e3);
Finterp_displacement = map2griddedInterpolant(map_img);
Hmax = max(map_img(:));
Hmin = min(map_img(:));
rays_span_angle = acos((Rbody + Hmin)/(Rbody + Hmax));  % Maximum span angle
fprintf(['\nMaximum ray span angle: ', num2str(rad2deg(rays_span_angle)), ' deg'])

% Sample sphere sectors
if phase_angle >= 0
    lon_CSF = linspace(phase_angle-pi/2-rays_span_angle, pi/2+rays_span_angle, round(sqrt(np)));
else
    lon_CSF = linspace(-pi/2-rays_span_angle, phase_angle+pi/2+rays_span_angle, round(sqrt(np)));
end
lat_CSF = linspace(-pi/2, pi/2, round(sqrt(np)));
[lon_CSF_grid, lat_CSF_grid] = meshgrid(lon_CSF, lat_CSF);
lon_CSF_all = lon_CSF_grid(:)';
lat_CSF_all = lat_CSF_grid(:)';
fprintf(['\nMin/max longitudes in CSF: ', num2str(rad2deg(min(lon_CSF_all))), ' / ', num2str(rad2deg(max(lon_CSF_all))), ' deg'])

% Extract positions and directions
pos_body2sec_CSF_ideal_all = cart_coord([repmat(Rbody, 1, np); lon_CSF_all; lat_CSF_all]);
sph_body2sec_IAU_ideal_all = sph_coord(dcm_CSF2IAU*pos_body2sec_CSF_ideal_all);
lon_IAU_all = sph_body2sec_IAU_ideal_all(2,:);
lat_IAU_all = sph_body2sec_IAU_ideal_all(3,:);
Rbody_vec = Rbody + Finterp_displacement(lat_IAU_all, lon_IAU_all);
pos_body2sec_CSF_all = cart_coord([Rbody_vec; lon_CSF_all; lat_CSF_all]);
pos_star2sec_CSF_all = pos_body2sec_CSF_all - pos_body2star_CSF; 
pos_cam2sec_CSF_all = pos_body2sec_CSF_all - pos_body2cam_CSF;
dir_star2sec_CSF_all = pos_star2sec_CSF_all./vecnorm(pos_star2sec_CSF_all);
dir_cam2sec_CSF_all = pos_cam2sec_CSF_all./vecnorm(pos_cam2sec_CSF_all);

% Find points close to terminator
ixs_in_term = find_sphere_terminator(lon_CSF_all, lat_CSF_all, rays_span_angle, k_threshold, method_threshold);
nterm =  sum(ixs_in_term);
pos_body2sec_CSF_term = pos_body2sec_CSF_all(:, ixs_in_term);
dir_star2sec_CSF_term = dir_star2sec_CSF_all(:, ixs_in_term);
pos_body2sec_IAU_term = dcm_CSF2IAU*pos_body2sec_CSF_term;
dir_star2sec_IAU_term = dcm_CSF2IAU*dir_star2sec_CSF_term;
fprintf(['\nTerminator sectors: ', num2str(nterm), ' out of ', num2str(np), ''])

% Find points close to limb
ixs_in_limb = find_sphere_slice(lon_CSF_all, lat_CSF_all, phase_angle, rays_span_angle, k_threshold, method_threshold);
nlimb =  sum(ixs_in_limb);
pos_body2sec_CSF_limb = pos_body2sec_CSF_all(:, ixs_in_limb);
dir_cam2sec_CSF_limb = dir_cam2sec_CSF_all(:, ixs_in_limb);
pos_body2sec_IAU_limb = dcm_CSF2IAU*pos_body2sec_CSF_limb;
dir_cam2sec_IAU_limb = dcm_CSF2IAU*dir_cam2sec_CSF_limb;
fprintf(['\nLimb sectors: ', num2str(nlimb), ' out of ', num2str(np), ''])

%% Run function

%--TERMINATOR
tic
[ixs_occluded_term, Hray_term, Rray_term] = find_sphere_occlusions(Rbody, pos_body2sec_IAU_term, dir_star2sec_IAU_term, Finterp_displacement, ...
                                    nrs, rays_span_angle, nworkers, method_intersection, method_sampling);
fprintf(['\nFound ', num2str(sum(ixs_occluded_term)), 'points occluded from star out of ', num2str(nterm)])
toc

%--LIMB
tic
[ixs_occluded_limb, Hray_limb, Rray_limb] = find_sphere_occlusions(Rbody, pos_body2sec_IAU_limb, dir_cam2sec_IAU_limb, Finterp_displacement, ...
                                    nrs, rays_span_angle, nworkers, method_intersection, method_sampling);
fprintf(['\nFound ', num2str(sum(ixs_occluded_limb)), 'points occluded from camera out of ', num2str(nterm)])
toc

%% POSTPROCESS
% extract first occlusion height of terminator points
d_rays_term = Hray_term - Rray_term;
d_rays_occluded_term = d_rays_term(:, ixs_occluded_term);
noccl_term = length(d_rays_occluded_term);
ixs_first_occlusion_on_ray_term = nan(1, noccl_term);
d_first_occlusion_term = nan(1, noccl_term);
for ix = 1:noccl_term
    ixs_first_occlusion_on_ray_term(ix) = find(d_rays_occluded_term(:, ix)<0, 1);
    d_first_occlusion_term(ix) = d_rays_occluded_term(ixs_first_occlusion_on_ray_term(ix), ix);
end

% extract first occlusion height of terminator points
d_rays_limb = Hray_limb - Rray_limb;
d_rays_occluded_limb = d_rays_limb(:, ixs_occluded_limb);
noccl_limb = length(d_rays_occluded_limb);
ixs_first_occlusion_on_ray_limb = nan(1, noccl_limb);
d_first_occlusion_limb = nan(1, noccl_limb);
for ix = 1:noccl_limb
    ixs_first_occlusion_on_ray_limb(ix) = find(d_rays_occluded_limb(:, ix)<0, 1);
    d_first_occlusion_limb(ix) = d_rays_occluded_limb(ixs_first_occlusion_on_ray_limb(ix), ix);
end


%% PLOTS
% distribution of occlusions: index on ray and distances
figure(), grid on, hold on
histogram(ixs_first_occlusion_on_ray_term, 1:nrs,'EdgeColor','none')
histogram(ixs_first_occlusion_on_ray_limb, 1:nrs,'EdgeColor','none')
xlabel('Sample index along ray')
legend('Terminator','Limb')

figure(), grid on, hold on
histogram(d_first_occlusion_term,'EdgeColor','none')
histogram(d_first_occlusion_limb,'EdgeColor','none')
xline(Hmin-Hmax, 'r')
legend('Terminator','Limb','Min DEM - Max DEM')
xlabel('Radial intersecting distance [m]')

figure(), grid on, hold on
plot(ixs_first_occlusion_on_ray_term, d_first_occlusion_term,'o')
plot(ixs_first_occlusion_on_ray_limb, d_first_occlusion_limb,'o')
xlabel('Sample index along ray')
ylabel('Radial intersecting distance [m]')
legend('Terminator','Limb')

%% 3D - TERMINATOR
if plot_mask_only
    cdata = ~ixs_occluded_term;
else
    cdata = double(ixs_occluded_term);
    cdata(ixs_occluded_term) = -1e-3*d_first_occlusion_term;
    cdata(~ixs_occluded_term) = nan;
end

[X, Y, Z] = sphere(1e3);
figure()
grid on, hold on, camzoom(1.1), view([0, 0])
axis equal
surf(Rbody*X, Rbody*Y, Rbody*Z, 1, 'EdgeColor','none','FaceColor','black','FaceAlpha', 0.2)
scatter3(pos_body2sec_CSF_term(1,:), pos_body2sec_CSF_term(2,:), pos_body2sec_CSF_term(3,:), sm, cdata, 'filled')
ampl = 2*Rbody;
quiver3(0,0,0, ampl*dir_body2star_CSF(1),ampl*dir_body2star_CSF(2),ampl*dir_body2star_CSF(3),'LineWidth',2)
quiver3(0,0,0, ampl*dir_body2cam_CSF(1),ampl*dir_body2cam_CSF(2),ampl*dir_body2cam_CSF(3),'LineWidth',2)
colormap(cmap);
col = colorbar;
if plot_mask_only
    col.Ticks = [0 1];
    col.TickLabels = {'false','true'};
    col.Label.String = 'Visibility [y/n]';
else
    col.Label.String = 'Occluding distance [km]';
end
legend('Moon','Points at Terminator','Body2Sun','Body2Cam')
view([180, 0])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

% % compare with DEM map
% Rbody_vec_temp = Rbody_vec(ixs_in_term);
% dem_visual = Rbody_vec_temp - Rbody;
% figure()
% grid on, hold on, camzoom(1.1), view([0, 0])
% axis equal
% surf(Rbody*X, Rbody*Y, Rbody*Z, 1, 'EdgeColor','none','FaceColor','black','FaceAlpha', 0.2)
% scatter3(pos_body2sec_CSF_term(1,:), pos_body2sec_CSF_term(2,:), pos_body2sec_CSF_term(3,:), sm, 1e-3*dem_visual, 'filled')
% ampl = 2*Rbody;
% quiver3(0,0,0, ampl*dir_body2star_CSF(1),ampl*dir_body2star_CSF(2),ampl*dir_body2star_CSF(3),'LineWidth',2)
% quiver3(0,0,0, ampl*dir_body2cam_CSF(1),ampl*dir_body2cam_CSF(2),ampl*dir_body2cam_CSF(3),'LineWidth',2)
% colormap(cmap);
% col = colorbar;
% col.Label.String = 'Displacement [km]';
% %col.Ticks = [0 1];
% %col.TickLabels = {'negative','positive'};
% legend('Moon','Displacement','Body2Sun','Body2Cam')
% view([180, 0])
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')

%% 3D - LIMB
if plot_mask_only
    cdata = ~ixs_occluded_limb;
else
    cdata = double(ixs_occluded_limb);
    cdata(ixs_occluded_limb) = -1e-3*d_first_occlusion_limb;
    cdata(~ixs_occluded_limb) = nan;
end

[X, Y, Z] = sphere(1e3);
figure()
grid on, hold on, camzoom(1.1), view([0, 0])
axis equal
surf(Rbody*X, Rbody*Y, Rbody*Z, 1, 'EdgeColor','none','FaceColor','black','FaceAlpha', 0.2)
scatter3(pos_body2sec_CSF_limb(1,:), pos_body2sec_CSF_limb(2,:), pos_body2sec_CSF_limb(3,:), sm, cdata, 'filled')
ampl = 2*Rbody;
quiver3(0,0,0, ampl*dir_body2star_CSF(1),ampl*dir_body2star_CSF(2),ampl*dir_body2star_CSF(3),'LineWidth',2)
quiver3(0,0,0, ampl*dir_body2cam_CSF(1),ampl*dir_body2cam_CSF(2),ampl*dir_body2cam_CSF(3),'LineWidth',2)
colormap(cmap);
col = colorbar;
if plot_mask_only
    col.Ticks = [0 1];
    col.TickLabels = {'false','true'};
    col.Label.String = 'Visibility [y/n]';
else
    col.Label.String = 'Occluding distance [km]';
end
legend('Moon','Points at Limb','Body2Sun','Body2Cam')
view([180, 0])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

% % compare with DEM map
% Rbody_vec_temp = Rbody_vec(ixs_in_limb);
% dem_visual = Rbody_vec_temp - Rbody;
% figure()
% grid on, hold on, camzoom(1.1), view([0, 0])
% axis equal
% surf(Rbody*X, Rbody*Y, Rbody*Z, 1, 'EdgeColor','none','FaceColor','black','FaceAlpha', 0.2)
% scatter3(pos_body2sec_CSF_limb(1,:), pos_body2sec_CSF_limb(2,:), pos_body2sec_CSF_limb(3,:), sm, 1e-3*dem_visual, 'filled')
% ampl = 2*Rbody;
% quiver3(0,0,0, ampl*dir_body2star_CSF(1),ampl*dir_body2star_CSF(2),ampl*dir_body2star_CSF(3),'LineWidth',2)
% colormap(cmap);
% col = colorbar;
% col.Label.String = 'Displacement [km]';
% %col.Ticks = [0 1];
% %col.TickLabels = {'negative','positive'};
% legend('Moon','Displacement','Body2Sun')
% view([180, 0])
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')
% 
%% Run function iteratively to find optimal nrs

nrs_vec = [2 4 6 8 10 20 40 60 80 100];  % ray sampling points
ixs_occluded_term = zeros(length(nrs_vec), nterm);
for ix = 1:length(nrs_vec)
    tic
    ixs_occluded_term(ix, :) = find_sphere_occlusions(Rbody, pos_body2sec_IAU_term, dir_star2sec_IAU_term, Finterp_displacement, nrs_vec(ix), ...
                        rays_span_angle, nworkers, method_intersection, method_sampling);
    fprintf(['\nFound ', num2str(sum(ixs_occluded_term(ix,:))), ' out of ', num2str(nterm),' occluded sectors'])
    toc
end

%%
figure()
subplot(1,2,1)
grid on, hold on
plot(nrs_vec, 1e2*sum(ixs_occluded_term, 2)./nterm)
plot(nrs_vec, 1e2*sum(ixs_occluded_term, 2)./np)
xlabel('Number of samples in each ray')
ylabel('Occluded sectors [%]')
legend('Out of the high-longitude ones','Out of the total')

subplot(1,2,2)
grid on, hold on
plot(nrs_vec(2:end), 1e2*diff(sum(ixs_occluded_term, 2))./nterm)
plot(nrs_vec(2:end),  1e2*diff(sum(ixs_occluded_term, 2))./np)
xlabel('Number of samples in each ray')
ylabel('Absolute Improvement [%]')
legend('Out of the high-longitude ones','Out of the total')
