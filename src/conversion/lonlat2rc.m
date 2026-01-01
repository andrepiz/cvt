function [rows, cols] = lonlat2rc(lon, lat, lonMin, lonMax, latMin, latMax, u, v)
% Convert longitude latitude coordinates to row column coordinates of a
% texture map of given lonMin, lonMax, latMin, latMax domain and number of
% horizontal and vertical pixels u and v. Longitudes are assumed to be
% increasing from left to right and latitudes are assumed to be decreasing
% from up to bottom. This means that when:
%   lat = latMax, rows = 1; 
%   lat = latMin, rows = v;
%   lon = lonMax, cols = u;
%   lon = lonMin, cols = 1;

cols = round((lon - lonMin) / (lonMax - lonMin) * (u-1)) + 1;

rows = round((latMax - lat) / (latMax - latMin) * (v-1)) + 1;

end