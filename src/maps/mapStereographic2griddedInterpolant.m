function [F, hlon, hlat] = mapStereographic2griddedInterpolant(map, lon, lat, flag_plot)
% Create a gridded interpolant from a stereographic map.
% The map is assumed to enter in matrix format with increasing rows going
% from north to south (90 to -90 degrees) and increasing columns going 
% from west to east (-180 to 180 degrees).

% if ~exist('lon_lims','var')
%     lon_lims = [-pi, pi];
% end
% 
% if ~exist('lat_lims','var')
%     lat_lims = [-pi/2, pi/2];
% end
% 
% [v, u] = size(map, [1,2]);
% lonspan = lon_lims(2) - lon_lims(1);
% latspan = lat_lims(2) - lat_lims(1);
% 
% % Compute the latitude and longitude grid points of a map of size u, v
% hlon = lonspan/u;
% hlat = latspan/v;
% 
% lonspanvec = [lon_lims(1), linspace(lon_lims(1) + 1.5*hlon, lon_lims(2) - 1.5*hlon, u - 2), lon_lims(2)];
% latspanvec = [lat_lims(1), linspace(lat_lims(1) + 1.5*hlat, lat_lims(2) - 1.5*hlat, v - 2), lat_lims(2)];
% 
% F =  griddedInterpolant({latspanvec, lonspanvec}, flip(map, 1), 'linear','none');

% %-- Assume a domain for the map
% %lonOriginal = rend.body.maps.displacement.limits(1,:);
% %latOriginal = rend.body.maps.displacement.limits(2,:);
% %lonNew = lonOriginal + deg2rad(10);
% %latNew = latOriginal + deg2rad(5);
% lonNewDeg = [-180 180];
% latNewDeg = [-90 -89];
% lonNewRad = deg2rad(lonNewDeg);
% latNewRad = deg2rad(latNewDeg);

lon_flat = lon(:);
lat_flat = lat(:);

% Scattered interpolant (can it be removed?)
Fscat = scatteredInterpolant(lon_flat, lat_flat, map(:), ...
                             'natural', 'none');

% Interpolating vectors
lon_lims = [min(lon_flat), max(lon_flat)];
lat_lims = [min(lat_flat), max(lat_flat)];

% Singularity fixes
if lon_lims(1) < -deg2rad(179) 
    lon_lims(1) = -pi;
end
if lon_lims(2) > deg2rad(179)
    lon_lims(2) = pi;
end
if lat_lims(1) < -deg2rad(-89) 
    lat_lims(1) = -pi/2;
end
if lat_lims(2) > deg2rad(89)
    lat_lims(2) = pi/2;
end

lonspan = lon_lims(2) - lon_lims(1);
latspan = lat_lims(2) - lat_lims(1);

% Derive u and v depending on density of provided data
% v is defined as the number of times the difference between the minimum 
% latitude angles fits in the latitude boundary
[diffsLat] = minNeighborDiffs(abs(lat));
res_v = max(diffsLat);
v = ceil(latspan/res_v);
u = 2*sum(size(map));

% Longitude sampling vector is uniform
lon_vec = linspace(lon_lims(1), lon_lims(2), u);
hlon = lonspan/u;

% Latitude sampling vector has more samples at high lat
lat_min = min(lat_flat);
lat_max = max(lat_flat);
alpha = 3;                        % fixed bias (can parameterize if you like)
w = linspace(0,1,v).^alpha;
lat_vec = lat_min + (lat_max - lat_min) * w;
% 
% [lat1, lat2, hlat] = sample_sphere_projected_uniform(0, v - 2, lat_lims);
% lat_vec = [lat1(1), 0.5*lat1 + 0.5*lat2, lat2(end)];

[Lon_eq, Lat_eq] = meshgrid(lon_vec, lat_vec);

% Evaluate the interpolant on the non-uniform grid
Z_eq = Fscat(Lon_eq, Lat_eq);      % [v x u]

% hlon = lonspan/u;
% hlat = latspan/v;
% 
% lonspanvec = [lon_lims(1), ...
%     linspace(lon_lims(1) + 1.5*hlon, lon_lims(2) - 1.5*hlon, u - 2), ...
%     lon_lims(2)];
% 
% latspanvec = [lat_lims(1), ...
%     linspace(lat_lims(1) + 1.5*hlat, lat_lims(2) - 1.5*hlat, v - 2), ...
%     lat_lims(2)];

F = griddedInterpolant({lat_vec, lon_vec}, Z_eq, 'linear','none');

if flag_plot
    figure();
    imagesc(rad2deg(lon_vec), rad2deg(lat_vec), F.Values);
    axis xy; colorbar; axis equal
    xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
    title('Stereographic DEM on equirectangular grid (griddedInterpolant values)');
end

end

function [diffs, neighbors, iMin, jMin] = minNeighborDiffs(A)
%MINNEIGHBORDIFFS  Difference between matrix minimum and its neighbors.
%
% [diffs, neighbors, iMin, jMin] = minNeighborDiffs(A)
%
% Inputs:
%   A      m x n matrix (real numbers)
%
% Outputs:
%   diffs      Vector of differences: neighbor - min(A) for each neighbor
%   neighbors  Vector of neighbor values
%   iMin, jMin Row/col indices of the minimum element
%
% Example:
%   A = [3 1 4; 2 5 6; 1 8 7];
%   [d, n, i, j] = minNeighborDiffs(A);  % d = [2 3 5 4]

% Find global minimum index
[~, idx] = min(A(:));
[iMin, jMin] = ind2sub(size(A), idx);
minVal = A(iMin, jMin);

% Define neighbor window (3x3 clamped at edges)
iRange = max(1, iMin-1) : min(size(A,1), iMin+1);
jRange = max(1, jMin-1) : min(size(A,2), jMin+1);

% Extract neighborhood and compute differences
neighborhood = A(iRange, jRange);
neighbors = neighborhood(:);
diffs = neighbors - minVal;

% Optionally exclude the center (minimum itself)
center_i = iMin - iRange(1) + 1;
center_j = jMin - jRange(1) + 1;
center_idx = sub2ind(size(neighborhood), center_i, center_j);

% Remove center from results (difference would be zero anyway)
neighbors(center_idx) = [];
diffs(center_idx) = [];

end