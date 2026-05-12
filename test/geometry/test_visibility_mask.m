cvt_install()

flag_debug = true;
nthreads = 14;

%% inputs

nlon = 1200;
nlat = 600;

lon = linspace(-pi, pi, nlon); 
lat = linspace(-pi/2, +pi/2, nlat);


[longrid, latgrid] = meshgrid(lon, lat);

pos_c2t_TAR = [-30833357.2042197;
               -44770570.5694731;
               0];
%pos_c2t_TAR = [0;0;norm(pos_c2t_TAR)];
dcm_TAR2CAM = [-0.235424743507297	0.188290215331253	0.953478791963676;
               -0.783919512449713	0.543122151702287	-0.300812812446280;
               -0.574495662076112	-0.818269556595554	0.0197399833993635];
pos_s2t_TAR = [1e12;0;0];
%dcm_TAR2CAM = eye(3);
Rbody = 1737.4e3;

fov = deg2rad(4); % round fov
%fov = deg2rad([2, 4]); % square fov

%% run
%%---
tic
hideout_hw = get_sphere_visibility_masks(latgrid, longrid, pos_c2t_TAR, pos_s2t_TAR, dcm_TAR2CAM, Rbody, fov, 1, flag_debug);
t_elapsed_hw = toc
%%--

%%---PARALLEL COMPUTATION
p = gcp('nocreate');
if isempty(p)
    parpool
end
tic
hideout_parhw = get_sphere_visibility_masks(latgrid, longrid, pos_c2t_TAR, pos_s2t_TAR, dcm_TAR2CAM, Rbody, fov, nthreads, flag_debug);
t_elapsed_parhw = toc
%%--