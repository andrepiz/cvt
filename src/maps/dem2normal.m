function [normal, rgb] = dem2normal(filepath_dem, filepath_normal, Rbody, dem_scaling, frame_normal, nbit_normal, flag_remove_outlier, flag_parallel, flag_double_prec, flag_plot)

if ~exist('flag_remove_outlier','var')
    flag_remove_outlier = true;
end

dem_data = imread(filepath_dem);
if flag_double_prec
    dem_data = double(dem_data);
end
dem_map = dem_scaling*dem_data;

if flag_remove_outlier
    ixs_outlier = abs(dem_map) > 1e3*median(abs(dem_map(:)));
    dem_map(ixs_outlier) = 0;
    if sum(ixs_outlier(:)) > 0
        warning(['DEM set to 0 for ', num2str(sum(ixs_outlier(:))), ' outlier (absolute values larger than 100 times the median)'])
    end
end

lon = cast(linspace(-pi, pi, size(dem_map, 2)), class(dem_data)); % increasing columns are increasing longitudes
lat = cast(linspace(pi/2, -pi/2, size(dem_map, 1)), class(dem_data)); % increasing rows are decreasing latitudes
[longrid, latgrid] = meshgrid(lon, lat);

height = dem_map + find_triaxial_radius(longrid, latgrid, Rbody);

[normal] = height2normal_vec(longrid, latgrid, height, struct('frame',frame_normal));

rgb = map2rgb(normal, nbit_normal);

try
    tiffCorrect(filepath_dem, filepath_normal, rgb, nbit_normal);
catch
    warning("Couldn't save RGB image")
end

if flag_plot
    figure(), 
    hold on, axis equal, 
    imagesc(rad2deg(lon), rad2deg(lat), height - mean(Rbody))
    colorbar
    colormap('turbo')
    xlim([-180 180])
    ylim([-90 90])
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    title('Displacement Map')

    figure(), 
    hold on, axis equal, 
    imagesc(rad2deg(lon), rad2deg(lat), rgb)
    xlim([-180 180])
    ylim([-90 90])
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    title('Normal Map')
end