function [clusterLabels, nCluster] = dbscan_approx(rowcol, epsilon, minpts, res_px)

np = size(rowcol, 1);
rows = rowcol(:, 1);
cols = rowcol(:, 2); 

% Build a binary mask from your filtered points
mask = false(res_px(1), res_px(2));  % use your actual image size
idx = sub2ind(size(mask), rows, cols);
mask(idx) = true;

% "epsilon" equivalent: dilate to connect points within epsilon px of each other
se = strel('disk', ceil(epsilon));
maskDilated = imdilate(mask, se);

% Label connected components
CC = bwconncomp(maskDilated, 8);    % connectivity tree default is 8

% Filter out clusters smaller than minpts
clusterSizes = cellfun(@numel, CC.PixelIdxList);
validClusters = find(clusterSizes >= minpts);

% Map back to per-point cluster labels
clusterLabels = zeros(1, np);
labelMatrix = labelmatrix(CC);
for i = 1:np
    lbl = labelMatrix(rows(i), cols(i));
    if ismember(lbl, validClusters)
        clusterLabels(i) = lbl;
    end
end
nCluster = numel(validClusters);

end