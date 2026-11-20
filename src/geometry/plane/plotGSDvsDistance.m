function plotGSDvsDistance(h, fov_deg, nPixels, maxDistance, nPoints, gsdMax)
%PLOTGSDVSDISTANCE Plot ground sampling distance versus ground distance.
%
% INPUTS
%   h           Vector of camera heights above ground [m]
%   fov_deg     Horizontal FOV [deg]
%   nPixels     Horizontal image resolution [pixels]
%   maxDistance Maximum ground distance to plot [m]
%   nPoints     Number of points used for each curve
%   gsdMax      Optional maximum GSD [m/pixel].
%               If provided, a line is drawn connecting the points where
%               each GSD curve reaches gsdMax.
%
% Example:
%   plotGSDvsDistance([2 5 10 20], 90, 1280, 100, 1000, 1);
%
%   The last argument can be omitted:
%   plotGSDvsDistance([2 5 10 20], 90, 1280, 100, 1000);

%% Input checks

if nargin < 5
    nPoints = 1000;
end

if nargin < 6
    gsdMax = [];
end

h = h(:);

assert(all(h > 0), ...
    'Camera heights must be positive.');

assert(fov_deg > 0 && fov_deg < 180, ...
    'FOV must be between 0 and 180 degrees.');

assert(nPixels > 1, ...
    'Number of pixels must be larger than 1.');

assert(maxDistance > 0, ...
    'Maximum distance must be positive.');

if ~isempty(gsdMax)
    assert(isscalar(gsdMax) && gsdMax > 0, ...
        'gsdMax must be a positive scalar.');
end

%% Camera parameters

fov = deg2rad(fov_deg);

% Angular size of one pixel
pixel_angle = fov / nPixels;

%% Ground distance

x = linspace(0.01, maxDistance, nPoints);

%% Plot

figure;
hold on;
grid on;
box on;

%% Store GSD threshold points

x_gsdMax = nan(size(h));

%% Compute GSD for each height

for i = 1:length(h)

    % Look angle from horizontal
    alpha = atan2(h(i), x);

    % Angular boundaries of one pixel
    alpha_1 = alpha - pixel_angle/2;
    alpha_2 = alpha + pixel_angle/2;

    % Ground intersections of pixel boundaries
    x_1 = h(i) ./ tan(alpha_1);
    x_2 = h(i) ./ tan(alpha_2);

    % Exact pixel ground footprint
    GSD = abs(x_1 - x_2);

    % Valid points
    valid = isfinite(GSD) & ...
            GSD > 0 & ...
            x_1 > 0 & ...
            x_2 > 0;

    % Plot GSD curve
    plot(x(valid), GSD(valid), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('h = %.1f m', h(i)));

    %% Find distance corresponding to gsdMax

    if ~isempty(gsdMax)

        x_valid = x(valid);
        gsd_valid = GSD(valid);

        % Find first crossing
        idx = find(gsd_valid >= gsdMax, 1, 'first');

        if ~isempty(idx)

            if idx == 1

                x_gsdMax(i) = x_valid(idx);

            else

                % Linear interpolation between the two surrounding points
                x1 = x_valid(idx-1);
                x2 = x_valid(idx);

                g1 = gsd_valid(idx-1);
                g2 = gsd_valid(idx);

                x_gsdMax(i) = x1 + ...
                    (gsdMax-g1) * ...
                    (x2-x1)/(g2-g1);

            end

        end

    end

end

%% Plot GSD threshold

if ~isempty(gsdMax)

    % Horizontal GSD threshold
    yline(gsdMax, '--', ...
        'LineWidth', 1.2, ...
        'DisplayName', sprintf('GSD_{max} = %.3g m/pixel', gsdMax));

    %% Connect threshold points

    validThreshold = isfinite(x_gsdMax);

    if nnz(validThreshold) >= 2

        plot(x_gsdMax(validThreshold), ...
             repmat(gsdMax, nnz(validThreshold), 1), ...
             '-o', ...
             'LineWidth', 2, ...
             'MarkerSize', 6, ...
             'DisplayName', 'GSD limit locus');

    end

end

%% Plot formatting

xlabel('Distance from camera nadir [m]');
ylabel('Ground sampling distance [m/pixel]');
title(sprintf( ...
    'Ground Sampling Distance vs. Ground Distance\nFOV = %.1f°, Resolution = %d px', ...
    fov_deg, nPixels));
legend('Location', 'northwest');
xlim([0 maxDistance]);

end