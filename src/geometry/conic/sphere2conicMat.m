function conicMat = sphere2conicMat(pos_origin2sphere_REF, pos_sphere2cam_REF, dcm_REF2CAM, K, R)

% Moon center [xc; yc; zc] Reference frame
xc = pos_origin2sphere_REF(1);
yc = pos_origin2sphere_REF(2);
zc = pos_origin2sphere_REF(3);

% Sphere conics
A_sphere = [eye(3)/R^2, -pos_origin2sphere_REF/R^2; -pos_origin2sphere_REF'/R^2, ((xc^2 + yc^2 + zc^2)/R^2 -1)];

% Define camera pose
tSC = -dcm_REF2CAM*pos_sphere2cam_REF; % Vector from Camera to world (In Camera frame)
rotSC = dcm_REF2CAM;                    % Camera orientation (Rotation from World to Cam)

% Camera matrix
P = K*[rotSC, tSC];

% Compute circle on the image by projecting the sphere 
conicMat = inv(P*(A_sphere\P'));

end