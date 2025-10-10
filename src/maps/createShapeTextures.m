function [Rbody, filename_displacement, filename_normal] = createShapeTextures(shape, sizes, ppd, frame, flag_plot)

% function handles to remove pole caps over given latitudes
fun_latmin2displ = @(lat, lat_min, R) -R*(1 - cos(pi/2-abs(lat_min))./cos(pi/2-abs(lat))).*(lat <= lat_min & lat <= 0);
fun_latmax2displ = @(lat, lat_max, R) -R*(1 - cos(pi/2-lat_max)./cos(pi/2-lat)).*(lat >= lat_max & lat >= 0);

% function handles to slice the sphere between given longitude range up to
% a maximum latitude
fun_frontup2displ = @(lon, lat, lon_min, lon_max, lat_lim, R) -R.*(1 - cos(lat_lim)./cos(lat)).*((lon >= lon_min & lon <= lon_max) & lat <= lat_lim & lat >= 0);
fun_frontbottom2displ = @(lon, lat, lon_min, lon_max, lat_lim, R) -R.*(1 - cos(lat_lim)./cos(lat)).*((lon >= lon_min & lon <= lon_max) & lat >= lat_lim  & lat <= 0);
fun_rearup2displ = @(lon, lat, lon_min, lon_max, lat_lim, R) -R.*(1 - cos(lat_lim)./cos(lat)).*((lon <= lon_min | lon >= lon_max) & lat <= lat_lim & lat >= 0);
fun_rearbottom2displ = @(lon, lat, lon_min, lon_max, lat_lim, R) -R.*(1 - cos(lat_lim)./cos(lat)).*((lon <= lon_min | lon >= lon_max) & lat >= lat_lim & lat <= 0) ;

% function handles to create sides at given longitude and latitude limit
fun_front2displ = @(lon, lat, lon_lim, lat_lim, R) -R.*(1 - abs(cos(lat_lim).*cos(lon_lim)./(cos(lat).*cos(lon)))).*(lon >= -abs(lon_lim) & lon <= abs(lon_lim) & lat <= abs(lat_lim.*cos(lon)./cos(lon_lim)) & lat >= -abs(lat_lim.*cos(lon)./cos(lon_lim)));
fun_rear2displ = @(lon, lat, lon_lim, lat_lim, R) -R.*(1 - abs(cos(lat_lim).*cos(lon_lim)./(cos(lat).*cos(lon)))).*((lon <= -abs(lon_lim) | lon >= abs(lon_lim)) & lat <= abs(lat_lim.*cos(lon)./cos(lon_lim)) & lat >= -abs(lat_lim.*cos(lon)./cos(lon_lim)));
fun_left2displ = @(lon, lat, lon_lim, lat_lim, R) -R.*(1 - abs(cos(lat_lim).*sin(lon_lim)./(cos(lat).*sin(lon)))).*(lon >= -pi + abs(lon_lim) & lon <= -abs(lon_lim) & lat <= abs(lat_lim.*sin(lon)./sin(lon_lim)) & lat >= -abs(lat_lim.*sin(lon)./sin(lon_lim)));
fun_right2displ = @(lon, lat, lon_lim, lat_lim, R) -R.*(1 - abs(cos(lat_lim).*sin(lon_lim)./(cos(lat).*sin(lon)))).*(lon >= abs(lon_lim) & lon <= pi - abs(lon_lim) & lat <= abs(lat_lim.*sin(lon)./sin(lon_lim)) & lat >= -abs(lat_lim.*sin(lon)./sin(lon_lim)));

% Init grid
ngrid = ppd*180;
hlonlat = pi/ngrid;
lon_vec = -pi + hlonlat/2:hlonlat:pi-hlonlat/2;
lat_vec = -pi/2 + hlonlat/2:hlonlat:pi/2 - hlonlat/2;
[lonGrid, latGrid] = meshgrid(lon_vec, lat_vec);
        
dem = zeros(ngrid, 2*ngrid);

switch shape 

    case 'cylinder'
        
        % [a, b] = [Height, Diameter]
        Rbody = sqrt(sizes(1)^2 + sizes(2)^2)/2;
        szn = sizes/Rbody;
        lat_lim = acos(szn(2)/2);
        lon_lims = [-pi, pi];
        for ilon = 1:2*ngrid
            for ilat = 1:ngrid
                displ = fun_frontup2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), lon_lims(1), lon_lims(2), lat_lim, 1) + ...
                        fun_rearup2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), lon_lims(1), lon_lims(2), lat_lim, 1) + ...
                        fun_frontbottom2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), lon_lims(1), lon_lims(2), lat_lim, 1) + ...
                        fun_rearbottom2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), lon_lims(1), lon_lims(2), lat_lim, 1) + ...
                        fun_latmin2displ(latGrid(ilat, ilon), lat_lim, 1) + ...
                        fun_latmax2displ(latGrid(ilat, ilon), lat_lim, 1);
                dem(ngrid - ilat + 1, ilon) = displ;
            end
        end

    case {'cube','plate','stick'}

        % [a, b, c] = [front, side, height]
        % a = 2*R*cos(lat_lim)*cos(lon_lim);
        % b = 2*R*cos(lat_lim)*sin(lon_lim);
        % c = 2*R*sin(lat_lim);
        Rbody = sqrt(sum(sizes.^2/4));
        szn = sizes/Rbody;
        lat_lim = asin(szn(3)/2);
        lon_lim = acos(szn(1)/2/cos(lat_lim));

        dem = zeros(ngrid, 2*ngrid);
        for ilon = 1:2*ngrid
            for ilat = 1:ngrid
                displ_bottomup = fun_latmin2displ(latGrid(ilat, ilon), -lat_lim, 1) + ...
                                 fun_latmax2displ(latGrid(ilat, ilon), lat_lim, 1);
                displ_frontrear = fun_front2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), -lon_lim, lat_lim, 1) + ...
                                  fun_rear2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), lon_lim, lat_lim, 1);
                displ_leftright = fun_left2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), -lon_lim, lat_lim, 1) +...
                                  fun_right2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), lon_lim, lat_lim, 1); 
                dem(ngrid - ilat + 1, ilon) = min([displ_bottomup, displ_frontrear, displ_leftright]);
            end
        end

    case 'capsule'
        
        % [a, b, c] = [Radius, Lower Base Diameter, Upper Base Diameter]
        % a = R
        % b = 2*R*cos(lat_min)
        % c = 2*R*cos(lat_max)
        Rbody = sizes(1);
        szn = sizes/(2*Rbody);
        lat_min = acos(szn(2)/2);
        lat_max = acos(szn(3)/2);
        lon_lims = [-pi, pi];
        for ilon = 1:2*ngrid
            for ilat = 1:ngrid
                displ = fun_frontup2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), lon_lims(1), lon_lims(2), lat_max, 1) + ...
                        fun_rearup2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), lon_lims(1), lon_lims(2), lat_max, 1) + ...
                        fun_frontbottom2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), lon_lims(1), lon_lims(2), lat_min, 1) + ...
                        fun_rearbottom2displ(lonGrid(ilat, ilon), latGrid(ilat, ilon), lon_lims(1), lon_lims(2), lat_min, 1);
                dem(ngrid - ilat + 1, ilon) = displ;
            end
        end

    case 'hockey'

        % [a, b, c] = [Radius, Height of Floor, Height of Ceiling]
        % a = R
        % b = R*sin(lat_min)
        % c = R*sin(lat_max)
        Rbody = sizes(1);
        szn = sizes/Rbody;
        lat_min = asin(szn(2));
        lat_max = asin(szn(3));
        for ilon = 1:2*ngrid
            for ilat = 1:ngrid
                displ = fun_latmax2displ(latGrid(ilat, ilon), lat_max, 1) + fun_latmin2displ(latGrid(ilat, ilon), lat_min, 1);
                dem(ngrid - ilat + 1, ilon) = displ;
            end
        end


    otherwise
        error('Shape not supported!')
end

filename_displacement = [shape,num2str(Rbody),'_displacement_',num2str(ppd),'.tif'];
filename_normal = [shape,num2str(Rbody),'_normal_', frame,'_',num2str(ppd),'.tif'];

map2tif(filename_displacement, single(dem), 32);
dem2normal(filename_displacement, filename_normal, 1, 1, frame, 16, false, true, flag_plot, false);

end

