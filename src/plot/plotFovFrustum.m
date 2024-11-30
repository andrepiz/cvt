function p = plotFovFrustum(posCam_REF, dcm_CAM2REF, fov, h, col, alpha)
% h: height
% b: base

sh = 10;    %scale factor for fov length wrt height

% Apex point (top of the pyramid)
apx_CAM = [0, 0, 0];

% Base vertices
bx = sh*h*tan(fov(1)/2);
by = sh*h*tan(fov(2)/2);

bs_CAM = [-bx/2, -by/2, sh*h;
           bx/2, -by/2, sh*h;
           bx/2, by/2, sh*h;
          -bx/2, by/2, sh*h];

% Combine apex and base vertices
verts_CAM = [apx_CAM; bs_CAM];

% Define the faces of the pyramid using the vertices
fcs = [
    1, 2, 3;  % Side face 1
    1, 3, 4;  % Side face 2
    1, 4, 5;  % Side face 3
    1, 5, 2;  % Side face 4
    2, 3, 4;  % Base face 1
    4, 5, 2;  % Base face 2
];

verts_REF = (posCam_REF + dcm_CAM2REF*verts_CAM');

% PLOT
p = patch('Vertices', verts_REF', 'Faces', fcs, ...
          'FaceColor', col, 'FaceAlpha', alpha);

end