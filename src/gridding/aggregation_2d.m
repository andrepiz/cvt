function dBins = aggregation_2d(dCoordsRC, dVals, dLimsRC, chAggregation, dGranularity)
% AGGREGATION_2D Group scattered points defined with decimal
% coordinates and compute an aggregated value for each fixed-size 2D discrete 
% bin by using a specific aggregation function.
%
%   dBins = aggregation_2d(dCoordsRC, dVals, dLimsRC, ...)
%
% INPUTS:
%   dCoordsRC     - 2xN array of decimal rows and columns coordinates
%   dValues       - 1xN vector of values
%   dLimsRC       - 2x2 array of [row_min row_max; col_min col_max] limits
%                   (default: tight bounds on dCoordsRC)
%   chAggregation - Aggregation function: 'mean', 'median', 'min', 'max', or 'mode'
%   dGranularity  - Grid upsampling (default: 1)
%
% OUTPUTS:
%   dBins          - 2D grid of aggregated values (NaN where undefined)

arguments
    dCoordsRC           (:, :) double {ismatrix}
    dVals               (1, :) double {isvector}
    dLimsRC             (:, 2) double {ismatrix} = [floor(min(dCoordsRC,[],2)) + 1, ceil(max(dCoordsRC,[],2))];
    chAggregation       char                   = 'min'
    dGranularity        (1, 1) double {isscalar} = 1
end

% Preliminary checks
ui32CoordColSize = uint32(size(dCoordsRC, 2));

assert(dGranularity >= 1 && mod(dGranularity, 1) == 0, 'Please provide granularity as a scalar integer larger or equal to 1')
assert(size(dCoordsRC,1) == 2, 'Please provide coords as [2xN] vector')
assert(size(dVals,1) == 1 && size(dVals,2) == ui32CoordColSize, 'Please provide values as [1xN] vector, where N is the second dimension of coords')
assert(size(dLimsRC,1) == 2 && size(dLimsRC,2) == 2, 'Please provide limits as [row_min row_max; col_min col_max] array specifying the min and max bins for the coordinates of the corresponding dimension')

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

pts = dCoordsRC*dGranularity;
lims  = [(dLimsRC(:,1)-1)*dGranularity, dLimsRC(:, 2)*dGranularity];
nrows = lims(1,2) - lims(1,1);
ncols = lims(2,2) - lims(2,1);
dBins = zeros(nrows, ncols);
nbins = numel(dBins);

% Ignore out of limits points
ixs_inside = pts(1,:) >= lims(1,1) & pts(1,:) <= lims(1,2) & ...
             pts(2,:) >= lims(2,1) & pts(2,:) <= lims(2,2);
vals = dVals(ixs_inside);
nvals = length(vals);
if nvals == 0
    return
end

% pts_bounded
pts_bounded = pts(:, ixs_inside);
row_idx = ceil(pts_bounded(1, :));  % Point at [0.3, 0.7] will snap to row 1 and col 1
col_idx = ceil(pts_bounded(2, :));

% Linear indexing and aggregation
lin_idx = sub2ind([nrows, ncols], row_idx, col_idx);
binsVec = accumarray(lin_idx', vals', [nbins, 1], local_stat, NaN);
dBins = reshape(binsVec, [nrows, ncols]);

end
