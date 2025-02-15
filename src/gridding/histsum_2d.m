function bins = histsum_2d(dCoords, dValues, dLimits, dGranularity)
% HISTSUM_2D Group scattered points and sum their values into fixed-size 
% 2D discrete bins with a nearest-neighbour method. 
% The bins are created upsampling the limits with a specified granularity.

arguments
    dCoords double
    dValues double
    dLimits double = [floor(min(dCoords,[],2)), 1 + ceil(max(dCoords,[],2))]
    dGranularity (1,1) double = 1
end

% Preliminary checks
ui32CoordRowSize = uint32(size(dCoords, 1));
ui32CoordColSize = uint32(size(dCoords, 2));

assert(dGranularity >= 1, 'Please provide granularity as a scalar integer larger or equal to 1')
assert(size(dValues,1) == 1 && size(dValues,2) == ui32CoordColSize, 'Please provide values as [1xN] vector, where N is the second dimension of coords')
assert(size(dLimits,1) == ui32CoordRowSize && size(dLimits,2) == 2, 'Please provide limits as [Dx2] vector, where D is the first dimension of coords')

% Scale by granularity factor
coords_scaled = dCoords*dGranularity;
nx = dLimits(1,2)*dGranularity;
ny = dLimits(2,2)*dGranularity;
xmin = (dLimits(1,1)-1)*dGranularity;
ymin = (dLimits(2,1)-1)*dGranularity;

% Partitioning into same-size submatrices as histcounts2 does not allows grids
% larger than 1024 x 1024
xsubmax = 1024;
xsublims = xmin + round(linspace(0, nx/xsubmax, ceil(nx/xsubmax)+1)*xsubmax);
ysubmax = 1024;
ysublims = ymin + round(linspace(0, ny/ysubmax, ceil(ny/ysubmax)+1)*ysubmax);

bins = zeros(nx, ny);

% Binning
for ii = 1:length(xsublims) - 1
    for jj = 1:length(ysublims) - 1

        xsubedges = xsublims(ii):1:xsublims(ii+1);
        ysubedges = ysublims(jj):1:ysublims(jj+1);

        subixs = coords_scaled(1,:) > xsubedges(1) & coords_scaled(1,:) <= xsubedges(end) & ...
                 coords_scaled(2,:) > ysubedges(1) & coords_scaled(2,:) <= ysubedges(end);
        subcoords_scaled = coords_scaled(:, subixs);
        if isempty(subcoords_scaled)
            continue
        end

        [subcounts, ~, ~, subrows, subcols] = histcounts2(subcoords_scaled(1,:), subcoords_scaled(2,:), xsubedges, ysubedges);

        subbins = zeros(size(subcounts));

        subvalues = dValues(:, subixs);

        for ix = 1:length(subvalues)
            subbins(subrows(ix), subcols(ix)) = subbins(subrows(ix), subcols(ix)) + subvalues(ix);
        end
        subbinsixs = {xsubedges(2:end), ysubedges(2:end)};
        bins(subbinsixs{:}) = subbins;

    end

end


end

