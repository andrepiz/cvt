cvt_install()

%%
fov = deg2rad([160, 170]);
%fov = deg2rad(150)
nper = 100;

f = 0.004; 
muPixel = 18e-6;
res_px = 2*f/muPixel*tan(fov/2);
if length(res_px) == 1
    res_px = [res_px, res_px];
end
K = [f./muPixel 0                res_px(1)/2;...
     0             f./muPixel    res_px(2)/2;...
     0             0                1]; % [-] projection matrix

per = get_fov_perimeter(fov, nper, false);
nper = size(per, 2);

d_rnd = 5 + randn(1, nper);
points_along_per = per.*d_rnd;
points_proj = K*per;
uv_proj = points_proj([1:2],:)./points_proj(3, :);

% PLOT
figure(), grid on, hold on, axis equal
zvec = zeros(1, nper);
quiver3(zvec, zvec, zvec, per(1, :), per(2, :), per(3, :))
scatter3(points_along_per(1,:), points_along_per(2,:), points_along_per(3,:),'o')
xlabel('x'), ylabel('y'), zlabel('z'), view([10, 20])
title('FOV')

figure(), grid on, hold on, axis equal
xlabel('u'), ylabel('v')
xlim([0, res_px(1)])
ylim([0, res_px(2)])
title('FOV')
for ix = 1:nper
scatter(uv_proj(1,ix), uv_proj(2,ix),'o')
%pause(0.1)
end



