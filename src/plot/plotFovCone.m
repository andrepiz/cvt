function p = plotFovCone(posCam_CAM, dcm_REF2CAM, fov, h, b, col, alpha)
% h: height
% b: base

npc = 20 + 1;
sh = 1;    %scale factor for fov length wrt height
% create cylnder
B = tan(fov/2)*h + b; % baffle major base
SB = (sh*(B-b)+b)/B; %scale factor for major base given sh
r = linspace(b, SB*B, npc); %cone radius
[X,Y,Z] = cylinder(r);

% create patch
obj_REF = surf2patch(X,Y,Z*h*sh);
verts_REF = obj_REF.vertices';
verts_NEW = (posCam_CAM*ones(1, size(verts_REF,2)) + dcm_REF2CAM*verts_REF)';
obj_NEW = obj_REF;
obj_NEW.vertices = verts_NEW;

% PLOT
p = patch(obj_NEW, 'FaceAlpha', alpha, 'FaceColor', col);

end