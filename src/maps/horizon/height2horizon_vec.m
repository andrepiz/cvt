function horizon = height2horizon_vec(longrid, latgrid, height, params)
% height2horizon_vec
% Computes the local horizon elevation angle at each grid point using a 
% vectorized, tile-based approach.
%
% INPUTS:
%   longrid [v x u]        - Grid of longitudes (radians)
%   latgrid [v x u]        - Grid of latitudes (radians)
%   height  [v x u] or griddedInterpolant
%                          - Height map or interpolant (meters)
%   params.granularity     - Granularity factor for angular sampling
%   params.spanmax         - Maximum span angle for 
%
% OUTPUT:
%   horizon [v x u]        - Horizon map (maximum elevation angle at each point, in radians)

if ~isfield(params,'granularity')
    params.granularity = 1;
end
if ~isfield(params,'spanmax')
    params.spanmax = 'auto';
end

[res_v, res_u] = size(longrid);

% Get map bounds
lon_lims = [min(longrid(:)), max(longrid(:))];
lat_lims = [min(latgrid(:)), max(latgrid(:))];

% Resolution and interpolant setup
if isa(height, 'griddedInterpolant')
    Finterp_height = height;
    lonspan = lon_lims(2) - lon_lims(1);
    latspan = lat_lims(2) - lat_lims(1);
    hlon = lonspan / res_u;
    hlat = latspan / res_v;
else
    [Finterp_height, hlon, hlat] = map2griddedInterpolant(height, lon_lims, lat_lims);
end

if strcmp(params.spanmax,'auto')
    spanmax = acos(min(Finterp_height.Values(:)) / max(Finterp_height.Values(:)));
else
    spanmax = params.spanmax;
end

% Offset grid (local sampling pattern)
delta = params.granularity;
lon_offsets = -spanmax/2 : delta * hlon : spanmax/2;
lat_offsets = -spanmax/2 : delta * hlat : spanmax/2;
[lat_offsets_grid, lon_offsets_grid] = ndgrid(lat_offsets, lon_offsets);
nOffsets = numel(lat_offsets_grid);

% Flatten base coordinate grid
lonA = longrid(:);
latA = latgrid(:);
nA = numel(latA);
hA = Finterp_height(latA, lonA);

% Create tiles
tiling = nA*nOffsets / (10e3*10e3);
idxStart = 1:ceil(nA/tiling):nA;
idxStart(end+1) = nA + 1;  % to close final segment

horizon_flat = pi/2 * ones(nA, 1);  % initialize with horizon max value

% Loop over chunks
for ix = 1:length(idxStart) - 1
    idxsTemp = idxStart(ix):(idxStart(ix+1)-1);
    nTile = length(idxsTemp);

    % Base coordinates expanded to [nTile x nOffsets]
    latAred = latA(idxsTemp);
    lonAred = lonA(idxsTemp);

    latAred_exp = repmat(latAred, 1, nOffsets);
    lonAred_exp = repmat(lonAred, 1, nOffsets);

    lat_offsets_all = repmat(lat_offsets_grid(:)', nTile, 1);
    lon_offsets_all = repmat(lon_offsets_grid(:)', nTile, 1);

    % Apply offsets
    latB = latAred_exp + lat_offsets_all;
    lonB = lonAred_exp + lon_offsets_all;

    % Handle pole wrap-around
    overPole = latB > pi/2;
    latB(overPole) = pi - latB(overPole);
    lonB(overPole) = pi + lonB(overPole);
    lonB = wrapToPi(lonB);

    % Interpolate neighbor heights
    hB = Finterp_height(latB, lonB);

    % Elevation angle from each A to neighborhood
    elA = height2el(hA(idxsTemp), lonAred, latAred, hB, lonB, latB);

    % Assign maximum elevation angle (i.e., horizon)
    horizon_flat(idxsTemp) = max(elA, [], 2);

    % Display progress every chunk
    fprintf('\rProcessing tile %d of %d (%.1f%%)', ix, length(idxStart)-1, 100*ix/(length(idxStart)-1));
end

% Reshape to original size
horizon = reshape(horizon_flat, res_v, res_u);
end
