cvt_install()

%%
fov = deg2rad([15, 10]);
fov = deg2rad([70, 170]);
%fov = deg2rad([4, 10]);
%fov = deg2rad([90]);
Rbody = 1737.4;
pos_body2cam_TAR = 2*Rbody*[0.8; 0.7; 0];
pos_body2cam_TAR = 0.8*Rbody*[0.8; 0.7; 0];
posBody_TAR = Rbody*[1; 1; 1];
%posBody_TAR = [0;0;0];
posCam_TAR = posBody_TAR + pos_body2cam_TAR;
nper = 100;
nsph = 100;

[tngts, dirs] = get_sphere_tangency_circle(posCam_TAR, posBody_TAR, Rbody, nper);

%% Plot the sphere
k1 = Rbody;
k3 = Rbody/2;

[x,y,z] = sphere(nsph);
figure(); grid on, hold on, axis equal, view(-dir_cam2body_TAR), camzoom(2)
xlabel('x'), ylabel('y'), zlabel('z')
surf(k1*x + posBody_TAR(1), k1*y + posBody_TAR(2), k1*z + posBody_TAR(3), 'EdgeColor','k','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.5);
scatter3(tngts(1,:), tngts(2,:), tngts(3,:),'g','LineWidth',3)
quiver3(repmat(posCam_TAR(1,:), 1, nper), repmat(posCam_TAR(2,:), 1, nper), repmat(posCam_TAR(3,:), 1, nper), ...
       k3*dirs(1, :), k3*dirs(2, :), k3*dirs(3, :),'r','LineWidth',1,'MarkerSize',1)
