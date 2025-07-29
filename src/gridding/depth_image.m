function [valsPixel, maskValid] = depth_image(dCoords, dValues, dLimits, ...
                                               chAggregation, ...
                                               dGranularity, bAntialiasing, chFilter)
% DEPTH_IMAGE Rasterizes scattered depth values onto a regular image grid.
%
%   [valsPixel, maskValid] = depth_image(dCoords, dValues, dLimits, ...)
%
% INPUTS:
%   dCoords       - 2xN array of [row; col] coordinates (e.g., [y; x])
%   dValues       - 1xN vector of depth values
%   dLimits       - 2x2 array of [row_min row_max; col_min col_max] limits
%                   (default: tight bounds on dCoords)
%   chAggregation - Aggregation function: 'mean', 'median', 'min', 'max', or 'mode'
%   dGranularity  - Grid spacing/resolution (default: 1)
%   bAntialiasing - Apply antialiasing on downsampling (true/false)
%   chFilter      - Filter name for antialiasing (e.g., 'gaussian')
%
% OUTPUTS:
%   valsPixel   - 2D grid of aggregated depth values (NaN where undefined)
%   maskValid   - Logical mask of valid (non-NaN) pixels

arguments
    dCoords             (:, :) double {ismatrix}
    dValues             (1, :) double {isvector}
    dLimits             (:, 2) double {ismatrix} = ...
        [floor(min(dCoords,[],2)), 1 + ceil(max(dCoords,[],2))];
    chAggregation       char                   = 'min'
    dGranularity        (1, 1) double {isscalar} = 1
    bAntialiasing       (1, 1) logical          = false
    chFilter            char                   = 'gaussian'
end

% Select aggregation function
switch lower(chAggregation)
    case 'mean'
        local_stat = @mean;
    case 'median'
        local_stat = @median;
    case 'min'
        local_stat = @min;
    case 'max'
        local_stat = @max;
    case 'mode'
        local_stat = @mode;
    otherwise
        error("Unsupported local statistic: '%s'. Use 'mean', 'median', 'min', 'max', or 'mode'.", chAggregation);
end

lims  = [(dLimits(:,1)-1)*dGranularity, dLimits(:, 2)*dGranularity];

% Grid size from limits
row_min = lims(1,1); row_max = lims(1,2);
col_min = lims(2,1); col_max = lims(2,2);
nRows = round((row_max - row_min) / dGranularity);
nCols = round((col_max - col_min) / dGranularity);

% Convert coordinates to grid indices
r = dCoords(1,:)';  % row index
c = dCoords(2,:)';  % column index
col_idx = floor((c - col_min) / dGranularity) + 1;
row_idx = floor((r - row_min) / dGranularity) + 1;

% Remove points outside the bounds
validMask = row_idx >= 1 & row_idx <= nRows & ...
            col_idx >= 1 & col_idx <= nCols;
row_idx = row_idx(validMask);
col_idx = col_idx(validMask);
d = dValues(validMask);

% Linear indexing and aggregation
lin_idx = sub2ind([nRows, nCols], row_idx, col_idx);
valsPixelVec = accumarray(lin_idx, d, [nRows*nCols, 1], local_stat, NaN);
valsPixelFine = reshape(valsPixelVec, [nRows, nCols]);

% Optional downsampling if granularity > 1
if dGranularity == 1
    valsPixel = valsPixelFine;
else
    if bAntialiasing
        kern_antialiasing = gaussianKernel(dGranularity, sqrt(dGranularity), false);
    else
        kern_antialiasing = [];
    end
    valsPixel = downsamplingreconstruction(valsPixelFine, dGranularity, kern_antialiasing, chFilter, dGranularity);
    valsPixel(valsPixel < 0) = 0;  % Clamp negative values
end

% Validity mask
maskValid = ~isnan(valsPixel);

end
