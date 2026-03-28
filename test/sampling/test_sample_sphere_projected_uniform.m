%% Test case
cvt_install()

alpha_test = deg2rad(-20);
nhphi = 100;
philims = 0.97*pi/2*[-1 1];

[phi1_test, phi2_test, hphi_test, err_test] = sample_sphere_projected_uniform(alpha_test, nhphi, philims);

figure(), 
subplot(1,3,1)
grid on, hold on
if alpha_test >= 0
    yline(rad2deg(alpha_test - pi/2),'k--')
    yline(rad2deg(pi/2),'r--')
else
    yline(rad2deg(-pi/2),'k--')
    yline(rad2deg(alpha_test + pi/2),'r--')
end
yline(rad2deg(alpha_test),'g--')
plot(1:nhphi, rad2deg(phi1_test))
plot(1:nhphi, rad2deg(phi2_test))
xlabel('point [#]')
ylabel('longitude [deg]')
legend('Starting Longitude','Ending Longitude','Alpha Angle','Longitude Sampling phi1','Longitude Sampling phi2')

subplot(1,3,2)
grid on, hold on
plot( 1:nhphi, rad2deg(hphi_test))
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
