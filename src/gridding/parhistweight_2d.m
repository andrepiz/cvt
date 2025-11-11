function bins = parhistweight_2d(dCoords, dValues, dLimits, dGranularity, ...
                                i32Algorithm, dWindowSize, dGaussianSigma, ...
                                bSnapPointsInsideLimits, ui8Numthreads)%#codegen
% PARHISTWEIGHT_2D Parallelized version of HISTWEIGHT_2D
arguments
    dCoords        (2, :) double {ismatrix}
    dValues        (1, :) double {isvector}
    dLimits        (:, 2) double {isvector} = [floor(min(dCoords,[],2)), 1 + ceil(max(dCoords,[],2))];
    dGranularity   (1, 1) double {isscalar} = 1 % uint32?
    i32Algorithm   (1, :) int32           = 2 % area
    dWindowSize    (1, 1) double {isscalar} = 1;
    dGaussianSigma (1, 1) double {isscalar} = 1/3
    bSnapPointsInsideLimits (1, 1) logical = false;    
    ui8Numthreads  (1,1) uint8 {isscalar} = 4
end

% Preliminary checks
ui32CoordRowSize = uint32(size(dCoords, 1));
ui32CoordColSize = uint32(size(dCoords, 2));

assert(dGranularity >= 1 && mod(dGranularity, 1) == 0, 'Please provide granularity as a scalar integer larger or equal to 1')
assert(size(dValues,1) == 1 && size(dValues,2) == ui32CoordColSize, 'Please provide values as [1xN] vector, where N is the second dimension of coords')
assert(size(dLimits,1) == ui32CoordRowSize && size(dLimits,2) == 2, 'Please provide limits as [Dx2] vector, where D is the first dimension of coords')

% Define coordinate shifts of each neighbor
dNeighborShifts = permn(-dWindowSize:1:dWindowSize, ui32CoordRowSize);

% Creating pools
idx_Pools = ceil(linspace(0, double(ui32CoordColSize), double(ui8Numthreads + uint8(1))));
npools = length(idx_Pools) - 1;
coords_pools = cell(1, npools);
values_pools = cell(1, npools);

for idx = 1:npools
    idx_temp = idx_Pools(idx) + 1:idx_Pools(idx + 1);
    coords_pools{idx} = dCoords(:, idx_temp);
    values_pools{idx} = dValues(:, idx_temp);
end

parfor (idx = 1:npools, ui8Numthreads)
%for idx = 1:npools
    [idxs_all{idx}, wvals_all{idx}] = parhistweight_2d_singlecall(coords_pools{idx}, values_pools{idx}, dLimits, dGranularity, ...
                                                i32Algorithm, dNeighborShifts, dGaussianSigma,...
                                                bSnapPointsInsideLimits);
end

% Find sizes
nrows = dLimits(1,2)*dGranularity;
ncols = dLimits(2,2)*dGranularity;

% Concatenate idx and wvals
idxs_all = cat(1, idxs_all{:});
vals_all = cat(1, wvals_all{:});

% Sum the values belonging to same index
wvals = accumarray(idxs_all, vals_all, [nrows*ncols, 1]);
bins = reshape(wvals, nrows, ncols);

end

