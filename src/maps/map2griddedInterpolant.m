function F = map2griddedInterpolant(map, lon_lims, lat_lims)
% Create a gridded interpolant from a longitude/latitude map.
% The map is assumed to enter in matrix format with increasing rows going
% from north to south (90 to -90 degrees) and increasing columns going 
% from west to east (-180 to 180 degrees).

if ~exist('lon_lims','var')
    lon_lims = [-pi, pi];
end

if ~exist('lat_lims','var')
    lat_lims = [-pi/2, pi/2];
end

[u, v] = size(map, [1,2]);
lonspan = lon_lims(2) - lon_lims(1);
latspan = lat_lims(2) - lat_lims(1);

% Compute the latitude and longitude grid points of a map of size u, v
hlon = lonspan/v;
hlat = latspan/u;

% The map values are assumed to be centered on the middle of the pixel
% Note: spanning vectors for interpolants must be sorted
% Note: boundaries must still be at defined limits to avoid gaps
lonspanvec = [lon_lims(1), lon_lims(1) + 1.5*hlon:hlon:lon_lims(2) - 1.5*hlon, lon_lims(2)];
latspanvec = [lat_lims(1), lat_lims(1) + 1.5*hlat:hlat:lat_lims(2) - 1.5*hlat, lat_lims(2)];
[latgrid, longrid] =  ndgrid(latspanvec, lonspanvec);

% map must be flipped upside down so to have increasing altitude with
% increasing rows
map_mod = flip(map, 1);

F =  griddedInterpolant(latgrid, longrid, map_mod, 'linear','none');

% F(0, 0)
% F(deg2rad(0), deg2rad(90))
% F(deg2rad(89), deg2rad(0))

%map_mod = pagetranspose(map);
%F =  griddedInterpolant(longrid, latgrid, map_mod);

% lon0P = rad2deg(longrid(v/2, u/2));
% lat0P = rad2deg(latgrid(v/2, u/2));
% map0P = map_mod(v/2, u/2, :);
% F0P = reshape(F(longrid(v/2, u/2), latgrid(v/2, u/2)), 3, 1);
% lonEP = rad2deg(longrid(3*v/4, u/2));
% latEP = rad2deg(latgrid(3*v/4, u/2));
% mapEP = map_mod(3*v/4, u/2, :);
% FEP = reshape(F(longrid(3*v/4, u/2), latgrid(3*v/4, u/2)), 3, 1);
% lonNP = rad2deg(longrid(v/2, u));
% latNP = rad2deg(latgrid(v/2, u));
% mapNP = map_mod(v/2, u, :);
% FNP = reshape(F(longrid(v/2, u), latgrid(v/2, u)), 3, 1);
% F0P(1), map0P(1)
% FEP(2), mapEP(2)
% FNP(3), mapNP(3)


end