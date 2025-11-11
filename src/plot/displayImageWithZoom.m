function fh = displayImageWithZoom(img, locZoom, posWinNormalized, lw, fhname, fhpos, fh)
% DISPLAYIMAGEWITHZOOM Displays an image with a zoomed portion overlayed on the original image.
% 
% Parameters:
%   img - The image data.
%   locZoom  - A 4-element vector specifying the [x, y, width, height]
%             location of the zoomed portion in the data
%   posWinNormalized  - A 4-element vector specifying the [x, y, widht, height] 
%                       normalized position of the zoom box in the figure
% 

% Read and display the image
[imgHeight, imgWidth, ~] = size(img);
ar = imgHeight/imgWidth;

if ~exist('fh','var')
    if ~exist('fhpos','var')
        fhpos = [100,100,600,600*ar];
    end
    fh = figure('name',fhname, 'units','pixels','position', fhpos); 
end
axMain = gca();
set(axMain,'Position', [0, 0, 1, 1]);  % Fill entire figure
imshow(img);
colormap('gray')
hold on;

for ix = 1:size(locZoom, 1)
 
    locZoom_temp = locZoom(ix, :);
    posWinNormalized_temp = posWinNormalized(ix, :);

    % Validate locZoom input
    if size(locZoom_temp, 2) ~= 4
        error('locZoom must be a vector with four elements: [x, y, width, height].');
    end
    
    % Validate posWinNormalized input
    if size(posWinNormalized_temp, 2) ~= 4
        error('posWinNormalized must be a vector with four elements: [x, y, width, height].');
    end

    % Extract coordinates for the zoomed portion
    x = locZoom_temp(1);
    y = locZoom_temp(2);
    width = locZoom_temp(3);
    height = locZoom_temp(4);
    
    % Ensure the zoom region is within image boundaries
    if x < 1 || y < 1 || x + width > imgWidth || y + height > imgHeight
        error('zoomLoc dimensions exceed the image boundaries.');
    end
    
    % Link the zoomed box with red lines to the inset corners
    line(axMain, [x, posWinNormalized_temp(1) * imgWidth], [y, (1 - posWinNormalized_temp(2) - posWinNormalized_temp(4)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
    line(axMain, [x + width, (posWinNormalized_temp(1) + posWinNormalized_temp(3)) * imgWidth], [y, (1 - posWinNormalized_temp(2) - posWinNormalized_temp(4)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
    line(axMain, [x, posWinNormalized_temp(1) * imgWidth], [y + height, (1 - posWinNormalized_temp(2)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
    line(axMain, [x + width, (posWinNormalized_temp(1) + posWinNormalized_temp(3)) * imgWidth], [y + height, (1 - posWinNormalized_temp(2)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
    
    xr = [posWinNormalized_temp(1) * imgWidth; (posWinNormalized_temp(1) + posWinNormalized_temp(3)) * imgWidth ; posWinNormalized_temp(1) * imgWidth;  (posWinNormalized_temp(1) + posWinNormalized_temp(3)) * imgWidth];
    yr = [(1 - posWinNormalized_temp(2) - posWinNormalized_temp(4)) * imgHeight; (1 - posWinNormalized_temp(2) - posWinNormalized_temp(4)) * imgHeight; (1 - posWinNormalized_temp(2)) * imgHeight;  (1 - posWinNormalized_temp(2)) * imgHeight];
    
    % Draw rectangles to indicate the zoomed area
    rectangle(axMain, 'Position', [x, y, width, height], 'EdgeColor', 'r', 'LineWidth', lw);
    rectangle(axMain, 'Position', [min(xr), min(yr), max(diff(xr)), max(diff(yr))], 'EdgeColor', 'r', 'LineWidth', 2*lw);
    
    % Create an inset of the zoomed portion within the same figure
    zoomedPortion = imcrop(img, [x, y, width, height]);
    axes(fh, 'units','normalized','Position', posWinNormalized_temp);
    imshow(zoomedPortion);
    colormap('gray')
end

end