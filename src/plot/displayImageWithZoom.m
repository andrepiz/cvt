function fh = displayImageWithZoom(img, winLoc, winPos, lw, fhname, fhpos)
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
[imgHeight, imgWidth, ~] = size(img);
ar = imgHeight/imgWidth;
if ~exist('fhpos','var')
    fhpos = [100,100,600,600*ar];
end
fh = figure('name',fhname, 'units','pixels','position', fhpos); 
axMain = gca();
set(axMain,'Position', [0, 0, 1, 1]);  % Fill entire figure
imshow(img);
colormap('gray')
hold on;

for ix = 1:size(winLoc, 1)
 
    winLoc_temp = winLoc(ix, :);
    winPos_temp = winPos(ix, :);

    % Validate winLoc input
    if size(winLoc_temp, 2) ~= 4
        error('winLoc must be a vector with four elements: [x, y, width, height].');
    end
    
    % Extract coordinates for the zoomed portion
    x = winLoc_temp(1);
    y = winLoc_temp(2);
    width = winLoc_temp(3);
    height = winLoc_temp(4);
    
    % Ensure the zoom region is within image boundaries
    if x < 1 || y < 1 || x + width > imgWidth || y + height > imgHeight
        error('zoomLoc dimensions exceed the image boundaries.');
    end
    
    % Link the zoomed box with red lines to the inset corners
    line(axMain, [x, winPos_temp(1) * imgWidth], [y, (1 - winPos_temp(2) - winPos_temp(4)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
    line(axMain, [x + width, (winPos_temp(1) + winPos_temp(3)) * imgWidth], [y, (1 - winPos_temp(2) - winPos_temp(4)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
    line(axMain, [x, winPos_temp(1) * imgWidth], [y + height, (1 - winPos_temp(2)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
    line(axMain, [x + width, (winPos_temp(1) + winPos_temp(3)) * imgWidth], [y + height, (1 - winPos_temp(2)) * imgHeight], 'Color', 'r', 'LineWidth', lw);
    
    xr = [winPos_temp(1) * imgWidth; (winPos_temp(1) + winPos_temp(3)) * imgWidth ; winPos_temp(1) * imgWidth;  (winPos_temp(1) + winPos_temp(3)) * imgWidth];
    yr = [(1 - winPos_temp(2) - winPos_temp(4)) * imgHeight; (1 - winPos_temp(2) - winPos_temp(4)) * imgHeight; (1 - winPos_temp(2)) * imgHeight;  (1 - winPos_temp(2)) * imgHeight];
    
    % Draw rectangles to indicate the zoomed area
    rectangle(axMain, 'Position', [x, y, width, height], 'EdgeColor', 'r', 'LineWidth', lw);
    rectangle(axMain, 'Position', [min(xr), min(yr), max(diff(xr)), max(diff(yr))], 'EdgeColor', 'r', 'LineWidth', 2*lw);
    
    % Create an inset of the zoomed portion within the same figure
    zoomedPortion = imcrop(img, [x, y, width, height]);
    axes(fh, 'units','normalized','Position', winPos_temp);
    imshow(zoomedPortion);
    colormap('gray')
end

end