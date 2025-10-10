function [normal, rgb] = dem2normal(filepath_dem, filepath_normal, Rbody, dem_scaling, dem_limits, frame_normal, nbit_normal, flag_remove_outlier, flag_parallel, flag_plot, fitting_error_threshold)

if ~exist('dem_scaling','var')
    dem_scaling = 1;
end

if ~exist('dem_limits','var')
    dem_limits = [-pi, pi; -pi/2, pi/2];
end

if ~exist('flag_remove_outlier','var')
    flag_remove_outlier = true;
end

if ~exist('flag_parallel','var')
    flag_parallel = false;
end

if ~exist('flag_plot','var')
    flag_plot = false;
end

dem_data = dem_scaling*double(imread(filepath_dem));

if ~exist('fitting_error_threshold','var')
    md = mean(abs([reshape(diff(dem_data, [], 1), 1, []), reshape(diff(dem_data, [], 2), 1, [])]));
    fitting_error_threshold = 1e-5*md;
end

if flag_remove_outlier
    ixs_outlier = abs(dem_data) > 1e3*median(abs(dem_data(:)));
    dem_data(ixs_outlier) = 0;
    if sum(ixs_outlier(:)) > 0
        warning(['DEM set to 0 for ', num2str(sum(ixs_outlier(:))), ' outlier (absolute values larger than 100 times the median)'])
    end
    clear ixs_outlier
end

lon_lims = dem_limits(1, :);
lat_lims = dem_limits(2, :);

[v, u] = size(dem_data, [1, 2]);
lonspan = lon_lims(2) - lon_lims(1);
latspan = lat_lims(2) - lat_lims(1);

% Compute the latitude and longitude grid points of a map of size u, v
hlon = lonspan/u;
hlat = latspan/v;

lon = [lon_lims(1) + 0.5*hlon:hlon:lon_lims(2) - 0.5*hlon];
lat = [lat_lims(2) - 0.5*hlat:-hlat:lat_lims(1) + 0.5*hlat];

[longrid, latgrid] = meshgrid(lon, lat);

height = dem_data + find_triaxial_radius(longrid, latgrid, Rbody);
clear dem_data

max_memory_multiplier = numel(longrid) / (10e3*10e3);
if max_memory_multiplier > 1
    split_factor = ceil(max_memory_multiplier);
    [rowStart, rowEnd, colStart, colEnd] = submatrixLimits(height, split_factor);
    normal = zeros(v, u, 3,'like', height);
    col_margin = 1;
    row_margin = 1;
    c = 0;
    for ii = 1:split_factor
        for jj = 1:split_factor
            rowIdxs = max(1, rowStart(ii, jj) - row_margin):min(v, rowEnd(ii, jj) + row_margin);
            colIdxs = max(1, colStart(ii, jj) - col_margin):min(u, colEnd(ii, jj) + col_margin);
            if flag_parallel
                normal_temp = height2normal_vec(longrid(rowIdxs, colIdxs), latgrid(rowIdxs, colIdxs), height(rowIdxs, colIdxs), struct('frame',frame_normal,'fitting_error_threshold',fitting_error_threshold));
            else
                normal_temp = height2normal(longrid(rowIdxs, colIdxs), latgrid(rowIdxs, colIdxs), height(rowIdxs, colIdxs), struct('frame',frame_normal,'fitting_error_threshold',fitting_error_threshold));
            end
            subrowIdxs = 1+row_margin:length(rowIdxs)-row_margin;
            subcolIdxs = 1+col_margin:length(colIdxs)-col_margin;
            normal(rowIdxs(subrowIdxs), colIdxs(subcolIdxs), :) = normal_temp(subrowIdxs, subcolIdxs, :);
            c = c + 1;
            disp(['Computing tile ',num2str(c), ' of ', num2str(split_factor^2), ': rows ', num2str(rowIdxs(1)), '-', num2str(rowIdxs(end)), ', columns ', num2str(colIdxs(1)), '-', num2str(colIdxs(end))])
        end
    end
    % correct meridian boundary with the mean between the two closest meridians
    normal(:, [1, u], :) = repmat(0.5*normal(:, 2, :) + 0.5*normal(:, u-1, :), 1, 2, 1);   
    % correct pole boundary with the mean across the closest parallel
    normal(1, :, :) = repmat(mean(normal(2, :, :), 2), 1, u, 1);
    normal(v, :, :) = repmat(mean(normal(v-1, :, :), 2), 1, u, 1);
else
    if flag_parallel
        [normal] = height2normal_vec(longrid, latgrid, height, struct('frame',frame_normal,'fitting_error_threshold',fitting_error_threshold));
    else
        [normal] = height2normal(longrid, latgrid, height, struct('frame',frame_normal,'fitting_error_threshold',fitting_error_threshold));
    end
end

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
    xlim([rad2deg(dem_limits(1,:))])
    ylim([rad2deg(dem_limits(2,:))])
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    title('Displacement Map')

    figure(), 
    hold on, axis equal, 
    imagesc(rad2deg(lon), rad2deg(lat), rgb)
    xlim([rad2deg(dem_limits(1,:))])
    ylim([rad2deg(dem_limits(2,:))])
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    title('Normal Map')
end

end