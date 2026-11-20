
map = imread('ldem_4.tif');
F = map2griddedInterpolant(map, [-pi pi], [-pi/2, pi/2]);

camDir_IAU = vecnormalize([0.1; 0.1; 0.1]);
sunDir_IAU = vecnormalize([-0.2; -0.4; 0.6]);
sphCam_IAU = sph_coord(camDir_IAU);
sphSun_IAU = sph_coord(sunDir_IAU);

[q_IAU2CSF, dcm_IAU2CSF, ixs_flip] = csf(sunDir_IAU, camDir_IAU);

%
figure()
hold on
imagesc(rad2deg(F.GridVectors{2}), rad2deg(F.GridVectors{1}),F.Values)
scatter(rad2deg(sphCam_IAU(2)), rad2deg(sphCam_IAU(3)),100,'b+');
scatter(rad2deg(sphSun_IAU(2)), rad2deg(sphSun_IAU(3)),100,'yd');

dcm_CSF2IAU = dcm_IAU2CSF';
lonCSF = 0;
latCSF = 0;
V = rotateGriddedInterpolation(F, lonCSF, latCSF, dcm_CSF2IAU)

lonIAU = sphSun_IAU(2);
latIAU = sphSun_IAU(3);
V_check = rotateGriddedInterpolation(F, lonIAU, latIAU, eye(3))

