
rng(10)
cm_map = 'parula';

flag_plot = true;

R = 300;
v1 = 1;
v2 = 10;
xshift = 20;
yshift = -40;

n = 1e5;
% 2D scattered points
xcoord = xshift + R*rand(1, n);
ycoord = yshift + R*rand(1, n);
values = [v1*ones(1, floor(n/2)), v2*ones(1, n-floor(n/2))];

pointcloud = [xcoord; ycoord; values];
weighter=2;
dx = 5;

tic
[Zidw, Zmin, Zmax, densitymap] = points2grid(pointcloud', weighter, dx, flag_plot);
toc