function projectSpheresOnSensor(camPos, bodiesPos, qCam, radii, focalLength, resolution, varargin)
% projectSpheresOnSensor  Visualize the projection of spherical bodies onto the
%                         camera pixel array over time.
%
% SYNTAX:
%   projectSpheresOnSensor(camPos, bodiesPos, qCam, radii, focalLength, resolution)
%   projectSpheresOnSensor(camPos, bodiesPos, qCam, radii, focalLength, resolution, bgImage)
%
% INPUTS:
%   camPos      - [3 x N]       Camera position wrt inertial origin at each epoch
%   bodiesPos   - [3 x N x M]   Bodies positions wrt inertial origin at each epoch
%                               (M bodies, N time steps)
%   qCam        - [4 x N] or [4 x 1]
%                               Unit quaternion [qx;qy;qz;qw] rotating vectors FROM
%                               camera frame TO inertial frame at each epoch.
%                               A single [4x1] column is broadcast to all epochs.
%   radii       - [1 x M]       Radius of each body (same units as positions)
%   focalLength - [scalar]      Camera focal length [pixels]
%   resolution  - [1 x 2]       Sensor resolution [width_px, height_px]
%
% OPTIONAL INPUT:
%   bgImage     - [height_px x width_px] or [height_px x width_px x 3]
%                 Background image displayed behind the projections.
%                 Must match the sensor resolution exactly.
%                 Grayscale (2-D) or RGB (3-D, uint8/double) accepted.
%                 Pass [] or omit to use the default dark background.
%
% DESCRIPTION:
%   For each epoch t and body k the function:
%     1. Rotates the camera-to-body vector into camera frame using qCam.
%     2. Projects the sphere centre with a pinhole model.
%     3. Computes the projected ellipse instead of a circle, accounting for
%        off-axis foreshortening.
%     4. Draws the animated ellipse on the pixel array.
%
% EXAMPLE:
%   N = 60; t = linspace(0, 2*pi, N);
%   camPos = zeros(3, N);
%   M = 2;
%   bodiesPos = zeros(3, N, M);
%   bodiesPos(:,:,1) = [1e4*cos(t); 1e4*sin(t); 5e4*ones(1,N)];
%   bodiesPos(:,:,2) = [-2e4*ones(1,N); 1e3*sin(t); 8e4*ones(1,N)];
%   ang  = linspace(0, pi/6, N);
%   qCam = [zeros(1,N); zeros(1,N); sin(ang/2); cos(ang/2)];
%   projectSpheresOnSensor(camPos, bodiesPos, qCam, [1000,800], 2000, [1024,1024]);

%--------------------------------------------------------------------------
%% Optional argument: background image
%--------------------------------------------------------------------------
p = inputParser;
addOptional(p, 'bgImage', [], @(x) isempty(x) || isnumeric(x));
parse(p, varargin{:});
bgImage = p.Results.bgImage;

%--------------------------------------------------------------------------
%% Input validation
%--------------------------------------------------------------------------
[~, N] = size(camPos);
M = length(radii);

assert(size(bodiesPos,2) == N, ...
    'bodiesPos must have N columns (same number of epochs as camPos).');
assert(size(bodiesPos,3) == M, ...
    'bodiesPos third dimension must equal length(radii) (number of bodies).');
assert(numel(resolution) == 2, ...
    'resolution must be a 2-element vector [width_px, height_px].');
assert(focalLength > 0, 'focalLength must be positive.');
assert(size(qCam,1) == 4, 'qCam must have 4 rows [qx; qy; qz; qw].');

if size(qCam,2) == 1
    qCam = repmat(qCam, 1, N);
end
assert(size(qCam,2) == N, ...
    'qCam must have N columns, or be a single [4x1] constant quaternion.');

width_px  = resolution(1);
height_px = resolution(2);
cx = width_px  / 2;
cy = height_px / 2;

hasBackground = ~isempty(bgImage);
if hasBackground
    assert(size(bgImage,1) == height_px && size(bgImage,2) == width_px, ...
        'bgImage must have size [%d x %d] to match resolution.', height_px, width_px);
end

%--------------------------------------------------------------------------
%% Figure and axes setup
%--------------------------------------------------------------------------
colours = lines(M);

figure('Name','Sphere Projections on Pixel Array', ...
       'Color',[0.08 0.08 0.10], ...
       'Units','normalized', 'OuterPosition',[0.05 0.05 0.9 0.9]);

axBgColor = [0.04 0.04 0.07];
if hasBackground
    axBgColor = 'none';
end

ax = axes('Color',  axBgColor, ...
          'XColor', [0.70 0.70 0.70], 'YColor', [0.70 0.70 0.70], ...
          'GridColor', [0.28 0.28 0.28], 'GridAlpha', 0.45, ...
          'TickDir','out', 'FontSize',11, ...
          'DataAspectRatio',[1 1 1]);
hold(ax,'on');
xlim(ax, [0 width_px]);
ylim(ax, [0 height_px]);
set(ax,'YDir','reverse');

if hasBackground
    imagesc(ax, [0 width_px], [0 height_px], bgImage);
    xlim(ax, [0 width_px]);
    ylim(ax, [0 height_px]);
    set(ax,'YDir','reverse');
end

grid(ax,'on');  box(ax,'on');

xlabel(ax, 'u  [px]', 'Color',[0.85 0.85 0.85], 'FontSize',12);
ylabel(ax, 'v  [px]', 'Color',[0.85 0.85 0.85], 'FontSize',12);

rectangle('Parent',ax, 'Position',[0 0 width_px height_px], ...
          'EdgeColor',[0.50 0.50 0.50], 'LineWidth',1.8);

plot(ax, [cx cx], [0 height_px], '--', 'Color',[0.32 0.32 0.32], 'LineWidth',0.9);
plot(ax, [0 width_px], [cy cy],  '--', 'Color',[0.32 0.32 0.32], 'LineWidth',0.9);

hLeg = gobjects(1,M);
for k = 1:M
    hLeg(k) = plot(ax, NaN, NaN, '-', 'Color', colours(k,:), 'LineWidth',2, ...
                   'DisplayName', sprintf('Body %d  (R = %.0f)',k,radii(k)));
end
legend(ax, hLeg, 'Location','northeast', ...
       'TextColor',[0.85 0.85 0.85], ...
       'Color',[0.11 0.11 0.14], 'EdgeColor',[0.35 0.35 0.35]);

titleHandle = title(ax, 'Initialising...', 'Color',[0.92 0.92 0.92], 'FontSize',13);

%--------------------------------------------------------------------------
%% Pre-allocate animated graphics objects
%--------------------------------------------------------------------------
theta_ell = linspace(0, 2*pi, 181);

hPatch  = gobjects(1,M);
hCentre = gobjects(1,M);
hOOF    = gobjects(1,M);

for k = 1:M
    hPatch(k)  = fill(ax, NaN, NaN, colours(k,:), ...
        'FaceAlpha',0.18, 'EdgeColor',colours(k,:), 'LineWidth',1.6);
    hCentre(k) = plot(ax, NaN, NaN, '+', ...
        'Color',colours(k,:), 'MarkerSize',8, 'LineWidth',1.6);
    hOOF(k)    = text(ax, 10 + (k-1)*160, height_px-18, '', ...
        'Color',colours(k,:), 'FontSize',9, 'Visible','off');
end

%--------------------------------------------------------------------------
%% Animation loop
%--------------------------------------------------------------------------
for t = 1:N

    titleHandle.String = sprintf('Epoch  %d / %d', t, N);

    for k = 1:M

        r_body = bodiesPos(:, t, k);
        r_cam  = camPos(:, t);
        dr_cam = rotframe(r_body - r_cam, qCam(:, t));

        X_c = dr_cam(1);
        Y_c = dr_cam(2);
        Z_c = dr_cam(3);
        d   = norm(dr_cam);

        if Z_c <= 0
            set(hPatch(k),  'XData',NaN, 'YData',NaN);
            set(hCentre(k), 'XData',NaN, 'YData',NaN);
            hOOF(k).String  = sprintf('Body %d: behind camera', k);
            hOOF(k).Visible = 'on';
            continue
        end

        if d <= radii(k)
            set(hPatch(k),  'XData',NaN, 'YData',NaN);
            set(hCentre(k), 'XData',NaN, 'YData',NaN);
            hOOF(k).String  = sprintf('Body %d: inside sphere', k);
            hOOF(k).Visible = 'on';
            continue
        end

        [ell_u, ell_v, u_e, v_e] = projectSphereEllipse( ...
            X_c, Y_c, Z_c, d, radii(k), focalLength, cx, cy, theta_ell);

        set(hPatch(k),  'XData',ell_u, 'YData',ell_v);
        set(hCentre(k), 'XData',u_e,   'YData',v_e);

        u_min = min(ell_u);  u_max = max(ell_u);
        v_min = min(ell_v);  v_max = max(ell_v);

        in_fov = (u_max > 0) && (u_min < width_px) && ...
                 (v_max > 0) && (v_min < height_px);

        if ~in_fov
            hOOF(k).String  = sprintf('Body %d: out of FOV', k);
            hOOF(k).Visible = 'on';
        else
            hOOF(k).Visible = 'off';
        end

    end

    drawnow limitrate;
    pause(0.03);
end

titleHandle.String = sprintf('Epoch  %d / %d  [complete]', N, N);

end


function [ell_u, ell_v, u_e, v_e, a, b, phi] = projectSphereEllipse( ...
        X_c, Y_c, Z_c, d, R, f, cx, cy, theta)
% projectSphereEllipse  Exact projection of a 3-D sphere onto the focal plane.
%
% Computes the exact conic (ellipse in the usual case) obtained by projecting
% a sphere with a pinhole camera, using the dual quadric / dual conic method.
%
% INPUTS:
%   X_c,Y_c,Z_c - sphere centre in camera frame
%   d           - norm([X_c;Y_c;Z_c])
%   R           - sphere radius
%   f           - focal length [px]
%   cx,cy       - principal point [px]
%   theta       - parametric angles for sampling [rad]
%
% OUTPUTS:
%   ell_u,ell_v - perimeter coordinates [px]
%   u_e,v_e     - ellipse centre [px]
%   a,b         - semi-major and semi-minor axes [px]
%   phi         - ellipse rotation angle [rad]

% Dual quadric of the sphere in homogeneous camera coordinates
% p = [X_c; Y_c; Z_c];
% Q_star = [eye(3), p;
%           p.',    d^2 - R^2];

% Camera matrix
K = [f, 0, cx;
     0, f, cy;
     0, 0, 1];
% P = K * [eye(3), zeros(3,1)];
% 
% % Project to dual conic and convert to primal conic
% C_star = P * Q_star * P.';
% C = inv(C_star);
% C = C / C(3,3);
% C = (C + C.') / 2;

C = sphere2conicMat([0;0;0], -[X_c; Y_c; Z_c], eye(3), K, R);

% Ellipse centre
centre = -C(1:2,1:2) \ C(1:2,3);
u_e = centre(1);
v_e = centre(2);

% Shape matrix
val_c = [centre; 1].' * C * [centre; 1];
M = -C(1:2,1:2) / val_c;

% Eigen-decomposition -> axes and orientation
[V, D] = eig(M);
lam = diag(D);

[lam, idx] = sort(lam);
V = V(:, idx);

a = 1 / sqrt(lam(1));
b = 1 / sqrt(lam(2));
phi = atan2(V(2,1), V(1,1));

% Parametric ellipse
cp = cos(phi);
sp = sin(phi);

ell_u = u_e + a*cp.*cos(theta) - b*sp.*sin(theta);
ell_v = v_e + a*sp.*cos(theta) + b*cp.*sin(theta);

end