function blendedImage = overlayTextureMask(textureImage, lonMin, lonMax, latMin, latMax, lonVec, latVec)
% OVERLAYTEXTUREMASK Overlay red mask on texture pixels with zero polygon points
% Inputs:
%   textureImage - HxWx[3] image array where upper row is north, rightmost
%   column is east
%   minLon,maxLon - longitude bounds of texture (radians)
%   minLat,maxLat - latitude bounds of texture (radians)  
%   lonVec,latVec - vectors of polygon vertices (radians)
% Outputs:
%   blendedImage - masked image

[H, W, C] = size(textureImage);

% Wrap around lon/lat points
[lonVec_new, latVec_new] = lonLatWrapAround(lonVec, latVec);

% Exlude points outside texture domain
ixs_inside = lonVec_new > lonMin & lonVec_new < lonMax & ...
            latVec_new > latMin & latVec_new < latMax;

% Convert lat/lon to image row/col indices
[rows, cols] = lonlat2rc(lonVec_new(ixs_inside), latVec_new(ixs_inside), lonMin, lonMax, latMin, latMax, W, H);

% Start with full red overlay on entire image
if isa(textureImage, 'uint8')
    textureImage = double(textureImage) / 255;
end
if size(textureImage, 3) == 1
    textureImage = repmat(textureImage, [1 1 3]);
end

redOverlay = repmat(reshape([1 0 0], [1 1 3]), [H W 1]);
blendedImage = textureImage;

% Convert row/col indices to linear indices
linearIdx = sub2ind([H W], rows(:), cols(:));

% Create mask of pixels with polygon points (unique linear indices)
hasPoints = false(H, W);
hasPoints(linearIdx) = true;

% % Reveal original texture where polygon points exist (remove red overlay)
for c = 1:size(blendedImage, 3)
    % Where revealMask=true: show original texture (alpha=0 for red)
    % Where revealMask=false: keep red overlay (alpha=1 for red)  
    blendedImage(:,:,c) = textureImage(:,:,c) .* hasPoints + ...
                          redOverlay(:,:,c) .* (~hasPoints);
end

end
