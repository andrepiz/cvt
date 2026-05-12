cvt_install()

%% Inputs

phase_angle = deg2rad(20);
R = 3000e3;

np = 100;

%% Local elevation 
d = 10*R;
ang_tangency = find_sphere_tangent_angle(d, R);
lon = linspace(phase_angle-ang_tangency, phase_angle+ang_tangency, np);
lat = linspace(-ang_tangency, ang_tangency, np)';

elevation = elevationFromPointOnSphere(lon, lat, R, d, phase_angle);

figure(), 
surf(rad2deg(lon) , rad2deg(lat), rad2deg(elevation));
colorbar
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim(rad2deg([phase_angle-ang_tangency, phase_angle+ang_tangency]))
ylim(rad2deg([-ang_tangency, ang_tangency]))
view([0, 90])

%% Critical elevation
elevation_critical = 

