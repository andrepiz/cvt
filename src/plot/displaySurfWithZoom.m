function fh = displaySurfWithZoom(xx, yy, zz, winLoc, winPos, lw, fh)
% DISPLAYIMAGEWITHZOOM Displays an image with a zoomed portion overlayed on the original image.
% 
% Parameters:
%   imagePath - Path to the image file.
%   winLoc  - A 4-element vector specifying the [x, y, width, height]
%             location in the data of the zoomed portion.
%   winPos  - A 4-element vector specifying the [x, y, widht, height] 
%             position in the axes of the zoomed portion
% 
% Example:
%   displayImageWithZoom('example.jpg', [50, 50, 100, 100])
% 

% Read and display the image
figure(fh)
ax1 = gca();
set(ax1,'Position', [0.2, 0.2, 0.6, 0.6]);
surf(xx, yy, zeros(size(zz)), zz, 'EdgeColor','none')
view(0, 90)
set(gca,'YDir','reverse')
xlabel('u [px]')
ylabel('v [px]')

% Validate zoomRect input
if numel(winLoc) ~= 4
    error('winLoc must be a vector with four elements: [x, y, width, height].');
end

% Extract coordinates for the zoomed portion
x = winLoc(1);
y = winLoc(2);
w = winLoc(3);
h = winLoc(4);

% Ensure the zoom region is within image boundaries
[imgHeight, imgWidth, ~] = size(zz);
if x < 1 || y < 1 || x + w > imgWidth || y + h > imgHeight
    error('winLoc dimensions exceed the image boundaries.');
end

% Draw a rectangle to indicate the zoomed area
rectangle('Position', [x, y, w, h], 'EdgeColor', 'r', 'LineWidth', lw);

% Create an inset of the zoomed portion within the same figure
ax2 = axes(fh, 'units','normalized','Position', winPos);
surf(xx, yy, zeros(size(zz)), zz, 'EdgeColor','none')
view(0, 90)
set(gca,'YDir','reverse')
colormap(ax1.Colormap)
clim(clim(ax1))
xlim([x, x + w]);
ylim([y, y + h]);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);

% Convert ax2 position from figure coordinates to ax1 coordinates
pos1 = ax1.Position; % [x, y, width, height] in normalized figure coordinates
pos2 = ax2.Position;

% Compute ax2 position relative to ax1
winPosRel(1) = (pos2(1) - pos1(1)) / pos1(3);
winPosRel(2) = (pos2(2) - pos1(2)) / pos1(4);
winPosRel(3) = pos2(3) / pos1(3);
winPosRel(4) = pos2(4) / pos1(4);

xr(1) = winPosRel(1) * imgWidth;
xr(2) = (winPosRel(1) + winPosRel(3)) * imgWidth;
xr(3) = winPosRel(1) * imgWidth;
xr(4) = (winPosRel(1) + winPosRel(3)) * imgWidth;
yr(1) = (1 - winPosRel(2) - winPosRel(4)) * imgHeight;
yr(2) = (1 - winPosRel(2) - winPosRel(4)) * imgHeight;
yr(3) = (1 - winPosRel(2)) * imgHeight;
yr(4) = (1 - winPosRel(2)) * imgHeight;

line(ax1, [x, xr(1)], [y, yr(1)], 'Color', 'r', 'LineWidth', lw);
line(ax1, [x + w, xr(2)], [y, yr(2)], 'Color', 'r', 'LineWidth', lw);
line(ax1, [x, xr(3)], [y + h, yr(3)], 'Color', 'r', 'LineWidth', lw);
line(ax1, [x + w, xr(4)], [y + h, yr(4)], 'Color', 'r', 'LineWidth', lw);

rectangle(ax1, 'Position', [xr(1), yr(1), xr(2) - xr(1), yr(4) - yr(1)], 'EdgeColor', 'r', 'LineWidth', 2*lw);

end