function outImg = mosaicMap(imageFiles, shape)
% mosaic_by_shape: Create a mosaic with specified shape [rows, cols].
% Inputs:
%   imageFiles - cell array of image filenames
%   shape - [rows, cols] specifying the mosaic layout
% Output:
%   outImg - mosaic image

    numImages = numel(imageFiles);
    rows = shape(1);
    cols = shape(2);
    
    if numImages > rows * cols
        error('Number of images exceeds mosaic capacity');
    end
    
    % Read all images into a cell
    imgs = cell(numImages, 1);
    for i = 1:numImages
        imgs{i} = imread(imageFiles{i});
    end
    
    % Resize all images to the smallest width and height
    widths = cellfun(@(im) size(im,2), imgs);
    heights = cellfun(@(im) size(im,1), imgs);
    minWidth = min(widths);
    minHeight = min(heights);
    
    for i = 1:numImages
        if size(imgs{i},1) ~= minHeight || size(imgs{i},2) ~= minWidth
            imgs{i} = imresize(imgs{i}, [minHeight, minWidth]);
        end
    end

    % Fill empty slots with black images if needed
    totalSlots = rows * cols;
    blackImage = zeros(minHeight, minWidth, size(imgs{1},3), 'like', imgs{1});
    for i = numImages+1:totalSlots
        imgs{i} = blackImage;
    end
    
    % Arrange images in column-major order:
    % The index in imgs is mapped to (row, col) as:
    % row = mod(i-1, rows) + 1
    % col = floor((i-1)/rows) +1
    % We build cell array for each column
    
    mosaicCols = cell(cols,1);
    for c = 1:cols
        colImgs = imgs((c-1)*rows + (1:rows));
        mosaicCols{c} = cat(1, colImgs{:}); % concatenate vertically in each column
    end
    
    % Concatenate all columns horizontally
    outImg = cat(2, mosaicCols{:});
end
