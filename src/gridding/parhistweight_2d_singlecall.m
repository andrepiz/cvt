function [idxs_all, wvals_all] = parhistweight_2d_singlecall(dCoords, dValues, dLimits, dGranularity, ...
                                                    i32Algorithm, dNeighborShifts, dGaussianSigma, ...
                                                    bSnapPointsInsideLimits)%#codegen
% PARHISTWEIGHT_2D_SINGLECALL Single call function of parhistweight_2d
%
% INPUTS:
%   dCoords        (2, :) double {ismatrix}
%   dValues        (1, :) double {isvector}
%   dLimits        (:, 2) double {isvector}
%   dGranularity   (1, 1) double {isscalar}
%   i32Algorithm   (1, :) int32              Algorithm of gridding. Method ID: 0: Inv. squared, 1: Diff, 2: area, 3: Gaussian 
%   dNeighborShifts(:, 2) double {isscalar}  Coordinate shifts of each neighboring bins
%   dGaussianSigma (1, 1) double {isscalar}  Only for "gaussian" method
%
% OUTPUTS:
%   bins           (M1, M2)

arguments
    dCoords        (2, :) double {ismatrix}
    dValues        (1, :) double {isvector}
    dLimits        (:, 2) double {isvector} = [floor(min(dCoords,[],2)), 1 + ceil(max(dCoords,[],2))]
    dGranularity   (1, 1) double {isscalar} = 1
    i32Algorithm   (1, :) int32             = 2
    dNeighborShifts(:, 2) double {isscalar} = permn(-1:1, 2)
    dGaussianSigma (1, 1) double {isscalar} = 1/2
    bSnapPointsInsideLimits (1, 1) logical = false
end

% Preliminary checks
ui32CoordRowSize = uint32(size(dCoords, 1));
ui32CoordColSize = uint32(size(dCoords, 2));

assert(dGranularity >= 1 && mod(dGranularity, 1) == 0, 'Please provide granularity as a scalar integer larger or equal to 1')
assert(size(dValues,1) == 1 && size(dValues,2) == ui32CoordColSize, 'Please provide values as [1xN] vector, where N is the second dimension of coords')
assert(size(dLimits,1) == ui32CoordRowSize && size(dLimits,2) == 2, 'Please provide limits as [Dx2] vector, where D is the first dimension of coords')

% Scale points and limits
pts = dCoords*dGranularity;
lims  = [(dLimits(:,1)-1)*dGranularity, dLimits(:, 2)*dGranularity];
nrows = lims(1,2);
ncols = lims(2,2);
nneighbors = size(dNeighborShifts, 1); 

% Check points outside limits
if bSnapPointsInsideLimits
    % Enforce bounds for sector_idx
    pts_bounded = max(lims(:,1), min(lims(:,2), pts));
    vals = dValues;
else
    % Ignore out of limits points
    ixs_inside = pts(1,:) >= lims(1,1) & pts(1,:) <= lims(1,2) & ...
                 pts(2,:) >= lims(2,1) & pts(2,:) <= lims(2,2);
    pts_bounded = pts(:, ixs_inside);
    vals = dValues(ixs_inside);
end
npts = size(pts_bounded, 2);

%-- INDEXES

% Find centers of each bin 
centers = (0.5 * sign(pts_bounded - round(pts_bounded)) + round(pts_bounded)); 

% Each center is shifted up to its neighbor
shiftedcenters_vert = centers(1,:) + dNeighborShifts(:, 1); 
shiftedcenters_horz = centers(2,:) + dNeighborShifts(:, 2);

% Check points outside limits
if bSnapPointsInsideLimits
    % Enforce bounds for rows and cols
    shiftedcenters_vert = max(lims(1,1), min(lims(1,2), shiftedcenters_vert));
    shiftedcenters_horz = max(lims(2,1), min(lims(2,2), shiftedcenters_horz));
    mask_inside = true(size(shiftedcenters_vert));
else
    % Ignore out of limits points
    mask_inside = shiftedcenters_vert >= lims(1,1) & shiftedcenters_vert <= lims(1,2) & ...
                  shiftedcenters_horz >= lims(2,1) & shiftedcenters_horz <= lims(2,2);
end

% Each bin row and col is given by ceiling the coordinate
rows = ceil(shiftedcenters_vert); 
cols = ceil(shiftedcenters_horz);

% Find linear indexes of each point
idxs_all = sub2ind([nrows, ncols], rows(mask_inside), cols(mask_inside));

%-- WEIGHTS
% Find distance of each point from the center of its neighboring bins
drows_neighbors = pts_bounded(1,:) - centers(1,:) - dNeighborShifts(:, 1);
dcols_neighbors = pts_bounded(2,:) - centers(2,:) - dNeighborShifts(:, 2);

% Constants
dWindowSize = round(log(nneighbors) / log(3)) - 1;
dConstDiffMethod = sqrt(2) * (0.5 + abs(dWindowSize)); % Constant scaling required by diff method
dConstGaussianMethod1 = 1 / (2 * dGaussianSigma * dGaussianSigma) ;
%dConstGaussianMethod2 = (1 / pi) * dConstGaussianMethod1;

% Init weights and indexes
dw_mat = zeros(nneighbors, npts);

for ii = 1:nneighbors

    d_neighbors = [drows_neighbors(ii, :); dcols_neighbors(ii, :)];

    % Find weights for the current shift
    switch i32Algorithm
        case 0 % Inverse squared method
            d = vecnorm(d_neighbors, 2, 1);
            dw_mat(ii,:) = 1./(d.^2);

        case 1 % Diff method
            % 1 minus distance normalized over maximum distance
            d = vecnorm(d_neighbors, 2, 1);
            dw_mat(ii,:) = 1 - d./dConstDiffMethod;

        case 2 % Area method
            % Fraction of [1x1] box area going to each sector
            ixs_neighbor = all(abs(d_neighbors) < 1, 1);
            dw_mat(ii, ixs_neighbor) = (1 - abs(d_neighbors(1, ixs_neighbor))).*(1 - abs(d_neighbors(2, ixs_neighbor)));

        case 3 % Gaussian kernel
            % Apply gaussian PSF
            d2 = d_neighbors(1, :).^2 + d_neighbors(2, :).^2;
            %dw_mat(ii, :) = dConstGaussianMethod2 * exp(-dConstGaussianMethod1 * d2);
            dw_mat(ii, :) = exp(-dConstGaussianMethod1 * d2);
            
        otherwise
            assert(0)
    end
end

% Set to zero the weights associated to points outside
dw_mat(~mask_inside) = 0;

% Normalize the weights to their sum for energy conservation
ws = sum(dw_mat, 1);
dw_mat = dw_mat./ws;

% Create weighted values matrix
wvals_mat = dw_mat.*vals;

% Exclude the values associated to points outside
wvals_all = wvals_mat(mask_inside);

end

