function displayImageWithZoom(img, zoomLoc, boxLoc, lw, fhname, fhpos)
% DISPLAYIMAGEWITHZOOM Displays an image with a zoomed portion overlayed on the original image.
% 
% Parameters:
%   imagePath - Path to the image file.
%   zoomRect  - A 4-element vector specifying the [x, y, width, height]
%               of the zoomed portion.
% 
% Example:
%   displayImageWithZoom('example.jpg', [50, 50, 100, 100])
% 

% Read and display the image
[imgHeight, imgWidth, ~] = size(img);
ar = imgHeight/imgWidth;
if ~exist('fhpos','var')
    fhpos = [100,100,600,600*ar];
end
fh = figure('name',fhname, 'units','pixels','position', fhpos); 
set(gca,'Position', [0, 0, 1, 1]);  % Fill entire figure
imshow(img);
colormap('gray')
hold on;

% Validate zoomRect input
if numel(zoomLoc) ~= 4
    error('zoomRect must be a vector with four elements: [x, y, width, height].');
end

% Extract coordinates for the zoomed portion
x = zoomLoc(1);
y = zoomLoc(2);
width = zoomLoc(3);
height = zoomLoc(4);

% Ensure the zoom region is within image boundaries
if x < 1 || y < 1 || x + width > imgWidth || y + height > imgHeight
    error('zoomRect dimensions exceed the image boundaries.');
end

% Link the zoomed box with red lines to the inset corners
line([x, boxLoc(1) * imgWidth], [y, (1 - boxLoc(2) - boxLoc(4)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
line([x + width, (boxLoc(1) + boxLoc(3)) * imgWidth], [y, (1 - boxLoc(2) - boxLoc(4)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
line([x, boxLoc(1) * imgWidth], [y + height, (1 - boxLoc(2)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
line([x + width, (boxLoc(1) + boxLoc(3)) * imgWidth], [y + height, (1 - boxLoc(2)) * imgHeight], 'Color', 'r', 'LineWidth', lw);

xr = [boxLoc(1) * imgWidth; (boxLoc(1) + boxLoc(3)) * imgWidth ; boxLoc(1) * imgWidth;  (boxLoc(1) + boxLoc(3)) * imgWidth];
yr = [(1 - boxLoc(2) - boxLoc(4)) * imgHeight; (1 - boxLoc(2) - boxLoc(4)) * imgHeight; (1 - boxLoc(2)) * imgHeight;  (1 - boxLoc(2)) * imgHeight];

% Draw a rectangle to indicate the zoomed area
rectangle('Position', [x, y, width, height], 'EdgeColor', 'r', 'LineWidth', lw);

rectangle('Position', [min(xr), min(yr), max(diff(xr)), max(diff(yr))], 'EdgeColor', 'r', 'LineWidth', 2*lw);

% Create an inset of the zoomed portion within the same figure
zoomedPortion = imcrop(img, [x, y, width, height]);
axes(fh, 'units','normalized','Position', boxLoc);
imshow(zoomedPortion);
colormap('gray')


end