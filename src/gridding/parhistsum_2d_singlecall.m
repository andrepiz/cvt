function [idxs, vals] = parhistsum_2d_singlecall(dCoords, dValues, dLimits, dGranularity)
% PARHISTSUM_2D_SINGLECALL Single call function of parhistsum_2d

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
coords = dCoords*dGranularity;
nx = dLimits(1,2)*dGranularity;
ny = dLimits(2,2)*dGranularity;
xmin = (dLimits(1,1)-1)*dGranularity;
ymin = (dLimits(2,1)-1)*dGranularity;

% Partitioning into same-size submatrices smaller than 1024 px
% (histcounts2 does not allows grids larger than 1024 x 1024)
xsubmax = 1024;
xsublims = xmin + round(linspace(0, nx/xsubmax, ceil(nx/xsubmax)+1)*xsubmax);
ysubmax = 1024;
ysublims = ymin + round(linspace(0, ny/ysubmax, ceil(ny/ysubmax)+1)*ysubmax);

% Obtain limits of each submatrix
[xsublimsleftgrid, ysublimsleftgrid] = meshgrid(xsublims(1:end-1), ysublims(1:end-1));
[xsublimsrightgrid, ysublimsrightgrid] = meshgrid(xsublims(2:end), ysublims(2:end));

% Init 
subidxsunique = cell(1, numel(xsublimsrightgrid));
subvals = subidxsunique;
for ix = 1:numel(xsublimsleftgrid)

    % Find bin edges
    xsubedges = xsublimsleftgrid(ix):1:xsublimsrightgrid(ix);
    ysubedges = ysublimsleftgrid(ix):1:ysublimsrightgrid(ix);

    % Extract coordinates belonging to the bins
    subixs = coords(1,:) > xsubedges(1) & coords(1,:) <= xsubedges(end) & ...
             coords(2,:) > ysubedges(1) & coords(2,:) <= ysubedges(end);
    subcoords = coords(:, subixs);
    if isempty(subcoords)
        continue
    end

    % Grid the continuous coordinates to the subbins
    % to check if it can be substituted by ceil
    [~, ~, ~, subrows, subcols] = histcounts2(subcoords(1,:), subcoords(2,:), xsubedges, ysubedges);

    % Convert the rows and cols to linear indexes
    subidxs = sub2ind([nx, ny], subrows + xsublimsleftgrid(ix), subcols + ysublimsleftgrid(ix));
    
    % Find unique indexes and sum the values belonging to same index
    [subidxsunique{ix},  ~, subgroups] = unique(subidxs);
    subvals{ix} = splitapply(@sum, dValues(subixs), subgroups');
end

idxs = [subidxsunique{:}];
vals = [subvals{:}];

end
