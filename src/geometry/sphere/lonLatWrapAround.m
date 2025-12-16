function [lon_new, lat_new] = lonLatWrapAround(lon, lat)
% Given longitude and latitude coordinates, fix their values so they are
% located on a equirectangular projection domain.

lon = wrapToPi(lon); % first set domain [-pi pi]
lat = wrapToPi(lat); % first set domain [-pi pi]

lon_new = lon;
lat_new = lat;

% wrap around points northern that north pole
lat_new(lat > pi/2) = pi - lat(lat > pi/2);
lon_new(lat > pi/2) = pi + lon(lat > pi/2);

% wrap around points southern that south pole
lat_new(lat < -pi/2) = -(lat(lat < -pi/2) + pi);
lon_new(lat < -pi/2) = pi + lon(lat < -pi/2);

lon_new = wrapToPi(lon_new); % set again domain [-pi pi]

end