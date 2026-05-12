% R = 1.4;
% th_vec = linspace(0, 2*pi, 50);
% pts = R*[cos(th_vec); -sin(th_vec)];
% nan_rand = randn(1, size(pts, 2))>0.2;
% pts(:, nan_rand) = nan;
% 

pts = load("pts.mat",'P_inter_REF');

% RUN
pts_filled = fill_nans_with_arcs(pts);

% PLOT

figure()
grid on, hold on, axis equal
plot(pts(1,3), pts(2,3),'s','MarkerSize',50)
plot(pts(1,end-2), pts(2,end-2),'s','MarkerSize',50)
plot(pts(1,:), pts(2,:),'ob-','LineWidth',2)
plot(pts_filled(1,:), pts_filled(2,:),'*k--','LineWidth',1)
legend('Start','End','Original','Filled')
