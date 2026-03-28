%% Test case
cvt_install()

fov = deg2rad(10);
%fov = deg2rad([10 15]);
f = 50e-3;
muPixel = 18e-6;
res_px = f/muPixel*2*tan(fov/2);
Rbody = 1737.4e3;
h = 1000e3;

distance = Rbody + h;
radius = Rbody;
offpoint = deg2rad(0);
phase_angle = deg2rad(30);
nhphi = 1e2;
rays_span_angle = deg2rad(8);

[ang_tangency, ang_bearing] = find_sphere_tangent_angle(distance, radius);
ang_tangency_modified = min(ang_tangency + rays_span_angle, pi/2);
phase_angle_modified = max(0, phase_angle - rays_span_angle);

philims = [-ang_tangency_modified, ang_tangency_modified];


[phi1, phi2, hphi] = sample_sphere_ifov_uniform(phase_angle_modified, distance, radius, offpoint, nhphi, philims, true);

%%
figure(), 
subplot(1,3,1)
grid on, hold on
if phase_angle >= 0
    yline(rad2deg(phase_angle - pi/2),'k--')
    yline(rad2deg(pi/2),'r--')
else
    yline(rad2deg(-pi/2),'k--')
    yline(rad2deg(phase_angle + pi/2),'r--')
end
yline(rad2deg(phase_angle),'g--')
plot(1:nhphi, rad2deg(phi1))
plot(1:nhphi, rad2deg(phi2))
xlabel('point [#]')
ylabel('longitude [deg]')
legend('Starting Longitude','Ending Longitude','Alpha Angle','Longitude Sampling phi1','Longitude Sampling phi2')

subplot(1,3,2)
grid on, hold on
plot( 1:nhphi, rad2deg(hphi))
xlabel('interval [#]')
ylabel('longitude [deg]')
legend('Longitude Step')

%% Span phase angle

alpha_vec = -pi:pi/100:pi;
nalpha = length(alpha_vec);
nhphi = 1000;

phi1 = zeros(nalpha, nhphi);
phi2 = zeros(nalpha, nhphi);
hphi = zeros(nalpha, nhphi);
err = zeros(nalpha, 1);
for ix = 1:nalpha
    [phi1(ix, :), phi2(ix, :), hphi(ix, :), err(ix)] = sample_sphere_projected_uniform(alpha_vec(ix), nhphi);
end

figure(), 
grid on, hold on
% yline(rad2deg(phi_sol0),'k--')
% yline(rad2deg(phi_soln),'r--')
% yline(rad2deg(alpha_test),'g--')
surf(1:nhphi, rad2deg(alpha_vec), rad2deg(phi1),'EdgeColor','none')
xlabel('point [#]')
ylabel('phase angle [deg]')
legend('Longitude Sampling phi1')
view(45, 10)

figure(), 
grid on, hold on
surf(1:nhphi, rad2deg(alpha_vec), rad2deg(hphi),'EdgeColor','none')
xlabel('point [#]')
ylabel('phase angle [deg]')
zlabel('longitude [deg]')
legend('Longitude Step')
set(gca,'ZScale','log')
view(45, 10)

figure(), 
grid on, hold on
plot(rad2deg(alpha_vec), rad2deg(err))
xlabel('phase angle [deg]')
ylabel('error [deg]')
