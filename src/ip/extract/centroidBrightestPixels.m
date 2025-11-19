function [centroidLoc, dnThreshold, dnThresholdBackground, maskNoise] = centroidBrightestPixels(img, fltSize, noiseThreshold, nfrac, method, flag_plot)
% CENTROIDBRIGHTESTPIXELS Find centroids of clusters formed by the brightest pixels in an image.
%
% [centroidLoc, dnThreshold, dnThresholdBackground, maskNoise] = centroidBrightestPixels(img, fltSize, nfrac, noiseThreshold, method)
%
% INPUTS:
%   img            - 2D grayscale image array.
%   fltSize        - Box filter size (set to 1 or [] to skip filtering).
%   nfrac          - Fraction (0 < nfrac <= 1, e.g., 0.01 for top 1%) or fixed count (N, e.g., 100 for top 100).
%   noiseThreshold - Numeric threshold, or 'dynamic' for mean plus 3*std, or omit to use default.
%   method         - 'moments' (default, intensity-weighted centroid) or 'median' (median location).
%
% OUTPUTS:
%   centroidLoc           - 2 x nClusters array: [column; row] centroids of each detected cluster.
%   dnThreshold           - Intensity threshold used to select bright pixels.
%   dnThresholdBackground - Threshold used to characterize background noise.
%   maskNoise             - Logical mask marking low-intensity (noise) pixels.
%
% NOTES:
%   - The function first filters and thresholds the image, then selects the brightest pixels.
%   - Clusters are found using DBSCAN on [row, col] of selected pixels (with optional inclusion of intensity, if desired).
%   - Centroids are computed via intensity-weighted moments or the median spatial location.
%   - Parameters should be tuned for the specifics of your data.
%
% EXAMPLE USAGE:
%   [centroidLoc, dnThreshold] = centroidBrightestPixels(img, 3, 0.01, 'dynamic', 'moments');

if ~exist('fltSize','var')
    imgData = double(img);
else
    imgData = double(imboxfilt(img, [fltSize fltSize]));
end

if ~exist('noiseThreshold','var')
    dnThreshold = mean(imgData(:)) + 3*std(imgData(:));
elseif strcmp(noiseThreshold,'dynamic')
    dnThreshold = mean(imgData(:)) + 3*std(imgData(:));
elseif isnumeric(noiseThreshold)
    dnThreshold = noiseThreshold;
end

if ~exist('nfrac','var')
    nfrac = 1;
end

if ~exist('method','var')
    method = 'moments';
end

if ~exist('flag_plot','var')
    flag_plot = false;
end

% Thresholding
maskNoise = imgData < dnThreshold;
dnThresholdBackground = mean(imgData(maskNoise)) + 3*std(imgData(maskNoise));
imgData(maskNoise) = 0;

[sortVals, sortInds] = sort(imgData(:), 'descend');

% Selection
if isnumeric(nfrac)
    if nfrac <= 1
        N = round(nfrac * numel(img));
    else
        N = min(nfrac, numel(img));
    end
else
    error('nfrac must be a number between 0 and 1 OR a number between 1 and the number of total elements in the img')
end

brightInds = sortInds(1:N);
[rows, cols] = ind2sub(size(imgData), brightInds);
brightVals = sortVals(1:N);
            
ixsFilt = brightVals > 0;

brightValsFilt = brightVals(ixsFilt);
rowsFilt = rows(ixsFilt);
colsFilt = cols(ixsFilt);

% Apply DBSCAN clustering on the filtered spatial coordinates
% you could use intensity as "dimension" to group based also on this
% information.
epsilon = 2.5;       % Distance threshold (tune as needed)
minpts = 5;          % Minimum points for a cluster (tune as needed)
clusterLabels = dbscan([rowsFilt, colsFilt], epsilon, minpts);
nCluster = max(clusterLabels(:));

centroidLoc = nan(2, nCluster);

for ix = 1:nCluster

    ixsCluster = clusterLabels == ix;
    rowsCluster = rowsFilt(ixsCluster);
    colsCluster = colsFilt(ixsCluster); 
    
    switch method
    
        case 'moments'
            % First moments of image weighted on element values
            valsCluster = brightValsFilt(ixsCluster);
            valsClusterSquared = valsCluster.*valsCluster;
            I00 = sum(valsClusterSquared);
            I10 = sum(rowsCluster.*valsClusterSquared);
            I01 = sum(colsCluster.*valsClusterSquared);
            centroidLoc(:, ix) = [I01/I00, I10/I00];
    
        case 'median'
            % Median of locations of values above noise threshold
            centroidLoc(:, ix) = [median(colsCluster), median(rowsCluster)];
    
        otherwise
            error('centroiding method not recognized')
    
    end
end

if flag_plot
    figure('name','centroiding','units','pixels','Position',[100 100 300 size(img,2)/size(img,1)*300])
    set(gca(), 'Position',[0 0 1 1])
    if max(img(:))>255
        imagesc(img)
        colormap('jet')
        col = 'white';
    else
        imshow(img);
        col = 'red';
    end
    hold on, axis equal
    scatter(centroidLoc(1, :), centroidLoc(2, :), 150, '+','LineWidth',0.2,'MarkerEdgeColor',col)
    scatter(centroidLoc(1, :), centroidLoc(2, :), 200, 'o','LineWidth',0.5,'MarkerEdgeColor',col)
    %ylim([0, size(img, 1)])
    %xlim([0, size(img, 2)])
end

end