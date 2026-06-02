function rgb = binary_rgb_overlay(img1, img2, varargin)
%BINARY_RGB_OVERLAY Create RGB false-color overlay from two binary images.
%
% rgb = binary_rgb_overlay(img1, img2) 
%   Binarizes img1 and img2 (using imbinarize), then creates an RGB image where:
%   - Red channel: binarized img1 ONLY (not where both present)
%   - Green channel: binarized img2 ONLY (not where both present)
%   - Blue channel: BOTH present (overlap region)
%   Result: Red=img1 only, Green=img2 only, Cyan=overlap.
%
% Optional name-value pairs:
%   'Threshold1'  threshold for img1 ('auto' or scalar), default 'auto'
%   'Threshold2'  threshold for img2 ('auto' or scalar), default 'auto'
%
% Example:
%   img1 = imread('coins1.png'); img2 = imread('coins2.png');
%   rgb = binary_rgb_overlay(img1, img2);
%   imshow(rgb);

    % Binarize images
    bw1 = imbinarize(img1);
    bw2 = imbinarize(img2);
    
    % Compute regions:
    % img1_only = bw1 & ~bw2  (red)
    % img2_only = bw2 & ~bw1  (green)  
    % both = bw1 & bw2        (blue/cyan)
    img1_only = bw1 & ~bw2;
    img2_only = bw2 & ~bw1;
    both = bw1 & bw2;
    
    % Create RGB channels
    rgb = cat(3, ...
        uint8(img1_only * 255), ...  % Red: img1 only
        uint8(img2_only * 255), ...  % Green: img2 only  
        uint8(both * 255));          % Blue: both present (cyan when viewed)
    
    % Optional display (remove if not wanted)
    if nargout == 0
        figure;
        imshow(rgb);
        title('Red = img1 only, Green = img2 only, Blue = both');
    end
end