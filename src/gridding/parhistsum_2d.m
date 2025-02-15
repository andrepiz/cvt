function bins = parhistsum_2d(dCoords, dValues, dLimits, dGranularity, ui8Numthreads)
% PARHISTSUM_2D Parallelized version of HISTSUM_2D

arguments
    dCoords double
    dValues double
    dLimits double = [floor(min(dCoords,[],2)), 1 + ceil(max(dCoords,[],2))]
    dGranularity (1,1) double = 1
    ui8Numthreads  (1,1) uint8 {isscalar} = 4
end

% Preliminary checks
ui32CoordRowSize = uint32(size(dCoords, 1));
ui32CoordColSize = uint32(size(dCoords, 2));

assert(dGranularity >= 1, 'Please provide granularity as a scalar integer larger or equal to 1')
assert(size(dValues,1) == 1 && size(dValues,2) == ui32CoordColSize, 'Please provide values as [1xN] vector, where N is the second dimension of coords')
assert(size(dLimits,1) == ui32CoordRowSize && size(dLimits,2) == 2, 'Please provide limits as [Dx2] vector, where D is the first dimension of coords')

idx_Pools = ceil(linspace(1, double(ui32CoordColSize), double(ui8Numthreads + uint8(1))));
ui8NumOfPools = length(idx_Pools) - 1;

% Creating pools
coords_pools = cell(1, ui8NumOfPools);
values_pools = cell(1, ui8NumOfPools);

for idx = 1:ui8NumOfPools
    idx_temp = idx_Pools(idx):idx_Pools(idx + 1);
    coords_pools{idx} = dCoords(:, idx_temp);
    values_pools{idx} = dValues(:, idx_temp);
end

parfor (idx = 1:ui8NumOfPools, ui8Numthreads)
%for idx = 1:ui8NumOfPools
    [idxs{idx}, vals{idx}] = parhistsum_2d_singlecall(coords_pools{idx}, values_pools{idx}, dLimits, dGranularity);
end

nx = dLimits(1,2)*dGranularity;
ny = dLimits(2,2)*dGranularity;
bins = zeros(nx, ny);
bins([idxs{:}]) = [vals{:}];

end