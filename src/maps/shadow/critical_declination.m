function [dec, debug] = critical_declination(lon, lat, ra, height, normal, params)
% Find the critical elevation at specified longitude and latitude along the
% specified right ascension direction, given the height and normal maps as
% matrices or gridded interpolants

% Critical declination is the minimum declination that is encountered scanning
% along the line defined by the right ascension with the theta parameter, 
% which is the polar coordinate in the longitude-latitude height map.

if ~exist('params.np','var')
    params.np = 1e2;
end
if ~exist('params.flag_debug','var')
    params.flag_debug = true;
end

% Extract height and normal at lon, lat
if isa(height, 'GriddedInterpolant')
    h = height(lon, lat);
    n = vecnormalize(squeeze(normal(lon, lat)));
else
    ix_lon = 1 + round(lon/2*pi*size(height,2));
    ix_lat = size(height,1) - (round(lat/pi*size(height,1)));
    h = height(ix_lat, ix_lon);
    ix_lon = 1 + round(lon/2*pi*size(normal,2));
    ix_lat = size(normal,1) - (round(lat/pi*size(normal,1)));
    n = vecnormalize(squeeze(normal(ix_lat, ix_lon, :)));
end

% Co-elevation is the angle between normal and z direction
coel = acos(n(3));

% Init 
th_vec = linspace(0, pi/2, params.np); % Theta can span up to 90 deg
lon_vec = lon + th_vec.*sin(ra);
lat_vec = lat + th_vec.*cos(ra);
dec = 0;
debug = struct();

% Cycle
for ix = 1:length(th_vec)
    
    th_temp = th_vec(ix);
    lon_temp = lon_vec(ix);
    lat_temp = lat_vec(ix);
    
    % Check bounds of lon, lat
    lon_temp = wrapTo2Pi(lon_temp);
    lat_temp = max(lat_temp, abs(lat_temp));
    lat_temp = min(lat_temp, pi - abs(lat_temp - pi));

    % Extract height and normal at lon_temp, lat_temp
    if isa(height, 'GriddedInterpolant')
        h_temp = height(lon_temp, lat_temp);
    else
        ix_lon = 1 + round(lon_temp/2*pi*size(height,2));
        ix_lat = size(height,1) - (round(lat_temp/pi*size(height,1)));
        h_temp = height(ix_lat, ix_lon);
    end

    % Compute critical declination
    phi = atan(((h/h_temp) + cos(th_temp))/sin(th_temp));
    if coel - phi < dec
        dec = coel - phi;
        if params.flag_debug
            debug.ix = ix;
        end
    end

    if params.flag_debug
        debug.th(ix) = th_temp;
        debug.lon(ix) = lon_temp;
        debug.lat(ix) = lat_temp;
        debug.h(ix) = h_temp;
        debug.dec(ix) = coel - phi;
    end

end