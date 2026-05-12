%%
cvt_install()

%% Inputs
AU = 149597870707;
Rbody = 1737.4e3;
d_body2cam = 80e6;
d_body2cam = 1837.4e3;
d_body2star = AU;
phase_angle_vec = deg2rad([-135 -45 45 135]);
% phase_angle_vec = deg2rad([-90 0 90 180]);
% phase_angle_vec = deg2rad([-90 0 90 180]) - deg2rad(5);
nc = length(phase_angle_vec);
np = 1e4;   % sectors

cam_tangency_angle = find_sphere_tangent_angle(d_body2cam, Rbody);
star_tangency_angle = find_sphere_tangent_angle(d_body2star, Rbody);

method_threshold = 'extended';
span_angle = deg2rad(8);
k_threshold = 0.1;

% 3d plots
cmap = 'parula';
sm = 5;

np = round(sqrt(np)).^2;

ixs_term = false(nc, np);
ixs_limb = false(nc, np);
ixs_before_term = false(nc, np);
ixs_before_limb = false(nc, np);
ixs_after_term = false(nc, np);
ixs_after_limb = false(nc, np);
dir_body2cam_CSF = zeros(3, nc);
for ix = 1:nc
% Prepare data
phase_angle = phase_angle_vec(ix);
pos_body2cam_CSF = d_body2cam*[cos(phase_angle);sin(phase_angle);0];
pos_body2star_CSF = d_body2star*[1;0;0];

dir_body2cam_CSF(:, ix) = pos_body2cam_CSF./vecnorm(pos_body2cam_CSF);
dir_body2star_CSF = pos_body2star_CSF./vecnorm(pos_body2star_CSF);

lon_CSF = linspace(-pi, pi, round(sqrt(np)));
lat_CSF = linspace(-pi/2, pi/2, round(sqrt(np)));
[lon_CSF_grid, lat_CSF_grid] = meshgrid(lon_CSF, lat_CSF);
lon_CSF_all = lon_CSF_grid(:)';
lat_CSF_all = lat_CSF_grid(:)';
fprintf(['\nMin/max longitudes in CSF: ', num2str(rad2deg(min(lon_CSF_all))), ' / ', num2str(rad2deg(max(lon_CSF_all))), ' deg'])

% Extract positions and directions
pos_body2sec_CSF_all = cart_coord([repmat(Rbody, 1, np); lon_CSF_all; lat_CSF_all]);

% RUN FUNCTIONS
[ixs_term(ix, :), ixs_before_term(ix, :), ixs_after_term(ix, :)] = find_sphere_terminator(lon_CSF_all, lat_CSF_all, span_angle, k_threshold, method_threshold, star_tangency_angle);
nterm = sum(ixs_term(ix,:));
fprintf(['\nTerminator sectors: ', num2str(nterm), ' out of ', num2str(np), ''])

[ixs_limb(ix, :), ixs_before_limb(ix, :), ixs_after_limb(ix, :)] = find_sphere_slice(lon_CSF_all, lat_CSF_all, phase_angle, span_angle, k_threshold, method_threshold, cam_tangency_angle);
nlimb = sum(ixs_limb(ix,:));
fprintf(['\nLimb sectors: ', num2str(nlimb), ' out of ', num2str(np), ''])

end

%% 3D
figure('units','normalized','Position',[0.1 0.1 0.8 0.7])

for ix = 1:nc
 
subplot(2,2,ix)
grid on, hold on
axis equal

cdata1 = double(ixs_term(ix,:));
cdata2 = 2*double(ixs_limb(ix,:));
cdata3 = 3*double(ixs_after_term(ix,:));
cdata4 = 4*double(ixs_after_limb(ix,:));
cdata0 = double(cdata1 > 0 | cdata2 > 0 | cdata3 > 0 | cdata4 > 0);
cdata0(cdata0 == 1) = nan;
cdata1(cdata1 == 0) = nan;
cdata2(cdata2 == 0) = nan;
cdata3(cdata3 == 0) = nan;
cdata4(cdata4 == 0) = nan;

[X, Y, Z] = sphere(1e3);
surf(Rbody*X, Rbody*Y, Rbody*Z, 1, 'EdgeColor','none','FaceColor','black','FaceAlpha', 0.2)
scatter3(pos_body2sec_CSF_all(1,:), pos_body2sec_CSF_all(2,:), pos_body2sec_CSF_all(3,:), sm, cdata1, 'filled')
scatter3(pos_body2sec_CSF_all(1,:), pos_body2sec_CSF_all(2,:), pos_body2sec_CSF_all(3,:), sm, cdata2, 'filled')
scatter3(pos_body2sec_CSF_all(1,:), pos_body2sec_CSF_all(2,:), pos_body2sec_CSF_all(3,:), sm, cdata3, 'filled')
scatter3(pos_body2sec_CSF_all(1,:), pos_body2sec_CSF_all(2,:), pos_body2sec_CSF_all(3,:), sm, cdata4, 'filled')
scatter3(pos_body2sec_CSF_all(1,:), pos_body2sec_CSF_all(2,:), pos_body2sec_CSF_all(3,:), sm, cdata0, 'filled')
ampl = 2*Rbody;
quiver3(0,0,0, ampl*dir_body2star_CSF(1),ampl*dir_body2star_CSF(2),ampl*dir_body2star_CSF(3),'LineWidth',2)
quiver3(0,0,0, ampl*dir_body2cam_CSF(1, ix),ampl*dir_body2cam_CSF(2, ix),ampl*dir_body2cam_CSF(3, ix),'LineWidth',2)
cm = colormap(cmap);
colormap(flip(cm));
col = colorbar;
col.Ticks = [0 1 2 3 4];
col.TickLabels = {'illuminated','terminator','limb','after terminator','after limb'};
col.Label.String = 'Region';
legend('Moon','Terminator','Limb','After Terminator','After Limb','Illuminated','Body2Sun','Body2Cam')
view([45 45])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

end