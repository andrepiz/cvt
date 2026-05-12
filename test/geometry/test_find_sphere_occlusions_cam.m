%%
cvt_install()

%% Inputs
AU = 149597870707;
Rbody = 1737.4e3;
phase_angle = -deg2rad(120);
np = 1e6;   % sectors
nrs = 50;  % ray sampling points
pos_body2cam_CSF = 80e6*[cos(phase_angle);sin(phase_angle);0];
%dcm_CSF2IAU = euler_to_dcm([pi/4; pi/4; pi/4]);
dcm_CSF2IAU = euler_to_dcm([0; 0; 0]);
filename = 'ldem_4.tif';
method_intersection = 'sampling';
method_sampling = 'linspace';
nworkers = 12;

method_threshold = 'poles';
k_threshold = 0.2;

% 3d plots
plot_mask_only = true;
cmap = 'parula';
sm = 2;

%% Prepare data
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
pos_star2sec_CSF_all = pos_body2sec_CSF_all - pos_body2cam_CSF; 
dir_star2sec_CSF_all = pos_star2sec_CSF_all./vecnorm(pos_star2sec_CSF_all);

% Extract only points around limb
switch method_threshold
    case 'simple'
        threshold_after_terminator = rays_span_angle.*cos(lat_CSF_all).^(k_threshold);
        ixs_in_terminator = abs(lon_CSF_all) - pi/2 < threshold_after_terminator;
    case 'poles'
        threshold_terminator = (pi/2 - rays_span_angle).*cos(lat_CSF_all).^(k_threshold);
        ixs_in_terminator = min(abs(lon_CSF_all), abs(abs(lon_CSF_all) - pi)) > threshold_terminator;
end
nhl =  sum(ixs_in_terminator);
pos_body2sec_CSF = pos_body2sec_CSF_all(:, ixs_in_terminator);
dir_star2sec_CSF = dir_star2sec_CSF_all(:, ixs_in_terminator);
pos_body2sec_IAU = dcm_CSF2IAU*pos_body2sec_CSF;
dir_star2sec_IAU = dcm_CSF2IAU*dir_star2sec_CSF;
fprintf(['\nTerminator sectors: ', num2str(nhl), ' out of ', num2str(np), ''])

%% Run function

tic
[ixs_occluded, Hray, Rray] = find_sphere_occlusions(Rbody, pos_body2sec_IAU, dir_star2sec_IAU, Finterp_displacement, ...
                                    nrs, rays_span_angle, nworkers, method_intersection, method_sampling);
toc
fprintf(['\nFound ', num2str(sum(ixs_occluded)), ' out of ', num2str(nhl),' occluded sectors'])

% extract first occlusion height
d_rays = Hray - Rray;
d_rays_occluded = d_rays(:, ixs_occluded);
noccl = length(d_rays_occluded);
ixs_first_occlusion_on_ray = nan(1, noccl);
d_first_occlusion = nan(1, noccl);
for ix = 1:noccl
    ixs_first_occlusion_on_ray(ix) = find(d_rays_occluded(:, ix)<0, 1);
    d_first_occlusion(ix) = d_rays_occluded(ixs_first_occlusion_on_ray(ix), ix);
end

%% PLOTS
% distribution of occlusions: index on ray and distances
figure(), grid on, hold on
histogram(ixs_first_occlusion_on_ray, 1:nrs)
xlabel('Sample index along ray')

figure(), grid on, hold on
histogram(d_first_occlusion)
xline(Hmin-Hmax, 'r')
legend('Hray-R@ray','Min DEM - Max DEM')
xlabel('Radial intersecting distance [m]')

figure(), grid on, hold on
plot(ixs_first_occlusion_on_ray, d_first_occlusion,'o')
xlabel('Sample index along ray')
ylabel('Radial intersecting distance [m]')

%% 3D
dir_sec2star_CSF = -dir_star2sec_CSF;
dir_body2star_CSF = pos_body2cam_CSF./vecnorm(pos_body2cam_CSF);

if plot_mask_only
    cdata = ixs_occluded;
else
    cdata = double(ixs_occluded);
    cdata(ixs_occluded) = -1e-3*d_first_occlusion;
    cdata(~ixs_occluded) = nan;
end

[X, Y, Z] = sphere(1e3);
figure()
grid on, hold on, camzoom(1.1), view([0, 0])
axis equal
surf(Rbody*X, Rbody*Y, Rbody*Z, 1, 'EdgeColor','none','FaceColor','black','FaceAlpha', 0.2)
scatter3(pos_body2sec_CSF(1,:), pos_body2sec_CSF(2,:), pos_body2sec_CSF(3,:), sm, cdata, 'filled')
ampl = 2*Rbody;
quiver3(0,0,0, ampl*dir_body2star_CSF(1),ampl*dir_body2star_CSF(2),ampl*dir_body2star_CSF(3),'LineWidth',2)
colormap(cmap);
col = colorbar;
if plot_mask_only
    col.Ticks = [0 1];
    col.TickLabels = {'false','true'};
    col.Label.String = 'Occluding [y/n]';
else
    col.Label.String = 'Occluding distance [km]';
end
legend('Moon','Occluded Sectors','Body2Sun')
view([180, 0])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

Rbody_vec_temp = Rbody_vec(ixs_in_terminator);
dem_visual = Rbody_vec_temp - Rbody;

% compare with DEM map
figure()
grid on, hold on, camzoom(1.1), view([0, 0])
axis equal
surf(Rbody*X, Rbody*Y, Rbody*Z, 1, 'EdgeColor','none','FaceColor','black','FaceAlpha', 0.2)
scatter3(pos_body2sec_CSF(1,:), pos_body2sec_CSF(2,:), pos_body2sec_CSF(3,:), sm, 1e-3*dem_visual, 'filled')
ampl = 2*Rbody;
quiver3(0,0,0, ampl*dir_body2star_CSF(1),ampl*dir_body2star_CSF(2),ampl*dir_body2star_CSF(3),'LineWidth',2)
colormap(cmap);
col = colorbar;
col.Label.String = 'Displacement [km]';
%col.Ticks = [0 1];
%col.TickLabels = {'negative','positive'};
legend('Moon','Displacement','Body2Sun')
view([180, 0])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

%% Run function iteratively to find optimal nrs

nrs_vec = [2 4 6 8 10 20 40 60 80 100];  % ray sampling points
ixs_occluded = zeros(length(nrs_vec), nhl);
for ix = 1:length(nrs_vec)
    tic
    ixs_occluded(ix, :) = find_sphere_occlusions(Rbody, pos_body2sec_IAU, dir_star2sec_IAU, Finterp_displacement, nrs_vec(ix), ...
                        rays_span_angle, nworkers, method_intersection, method_sampling);
    fprintf(['\nFound ', num2str(sum(ixs_occluded(ix,:))), ' out of ', num2str(nhl),' occluded sectors'])
    toc
end

%%
figure()
subplot(1,2,1)
grid on, hold on
plot(nrs_vec, 1e2*sum(ixs_occluded, 2)./nhl)
plot(nrs_vec, 1e2*sum(ixs_occluded, 2)./np)
xlabel('Number of samples in each ray')
ylabel('Occluded sectors [%]')
legend('Out of the high-longitude ones','Out of the total')

subplot(1,2,2)
grid on, hold on
plot(nrs_vec(2:end), 1e2*diff(sum(ixs_occluded, 2))./nhl)
plot(nrs_vec(2:end),  1e2*diff(sum(ixs_occluded, 2))./np)
xlabel('Number of samples in each ray')
ylabel('Absolute Improvement [%]')
legend('Out of the high-longitude ones','Out of the total')
