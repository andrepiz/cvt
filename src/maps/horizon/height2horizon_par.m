function horizon = height2horizon_par(longrid, latgrid, height, params)
% height2horizon_par
% Optimized and parallelized computation of horizon elevation angle at each grid point.
%
% INPUTS:
%   longrid [v x u]         - Grid of longitudes (radians)
%   latgrid [v x u]         - Grid of latitudes (radians)
%   height  [v x u] or griddedInterpolant
%                          - Height map or interpolant (meters)
%   params.granularity      - Sampling granularity (default: 1)
%   params.spanmax          - Max angular span (default: auto)
%
% OUTPUT:
%   horizon [v x u]         - Horizon map (max elevation angle at each point)

if ~isfield(params, 'granularity')
    params.granularity = 1;
end
if ~isfield(params, 'spanmax')
    params.spanmax = 'auto';
end

[res_v, res_u] = size(longrid);
lon_lims = [min(longrid(:)), max(longrid(:))];
lat_lims = [min(latgrid(:)), max(latgrid(:))];

% Interpolant setup
if isa(height, 'griddedInterpolant')
    Finterp_height = height;
    lonspan = lon_lims(2) - lon_lims(1);
    latspan = lat_lims(2) - lat_lims(1);
    hlon = lonspan / res_u;
    hlat = latspan / res_v;
    height_vals = Finterp_height.Values(:);
else
    [Finterp_height, hlon, hlat] = map2griddedInterpolant(height, lon_lims, lat_lims);
    height_vals = height(:);
end

% Compute max span
if strcmp(params.spanmax, 'auto')
    spanmax = acos(min(height_vals) / max(height_vals));
else
    spanmax = params.spanmax;
end

% Offset grid (sampling pattern)
delta = params.granularity;
lon_offsets = -spanmax/2 : delta * hlon : spanmax/2;
lat_offsets = -spanmax/2 : delta * hlat : spanmax/2;
[lat_offsets_grid, lon_offsets_grid] = ndgrid(lat_offsets, lon_offsets);
lat_offsets_vec = lat_offsets_grid(:);
lon_offsets_vec = lon_offsets_grid(:);

% Flatten base coordinate grid
lonA = longrid(:);
latA = latgrid(:);
nA = numel(latA);
hA = Finterp_height(latA, lonA);

% Preallocate output
horizon_flat = pi/2 * ones(nA, 1);

% Parallel loop over all points
parfor ix = 1:nA
    % Coordinates of point A
    lat0 = latA(ix);
    lon0 = lonA(ix);
    h0 = hA(ix);

    % Neighbor positions (nOffsets x 1)
    latB = lat0 + lat_offsets_vec;
    lonB = lon0 + lon_offsets_vec;

    % Handle polar wrap-around
    overPole = latB > pi/2;
    latB(overPole) = pi - latB(overPole);
    lonB(overPole) = pi + lonB(overPole);
    lonB = wrapToPi(lonB);

    % Interpolate neighbor heights
    hB = Finterp_height(latB, lonB);

    % Compute elevation angles from A to B points
    elA = height2el(h0, lon0, lat0, hB, lonB, latB);

    % Store max elevation (horizon)
    horizon_flat(ix) = max(elA);
end

% Reshape to original grid size
horizon = reshape(horizon_flat, res_v, res_u);

end
