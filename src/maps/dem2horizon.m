function [horizon] = dem2horizon(filepath_dem, filepath_horizon, Rbody, dem_scaling, dem_limits, nbit_map, flag_remove_outlier, flag_parallel, flag_plot, granularity, spanmax, split_factor)

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

if ~exist('granularity','var')
    granularity = 1;
end

if ~exist('spanmax','var')
    spanmax = 'auto';
end

if ~exist('split_factor','var')
    split_factor = 'auto';
end

dem_data = dem_scaling*double(imread(filepath_dem));

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

if strcmp(spanmax,'auto')
    spanmax_nominal = acos(min(height(:))./max(height(:)));
else
    spanmax_nominal = spanmax;
end

if strcmp(split_factor,'auto')
    max_number_elems =  spanmax_nominal/(granularity * hlon) * spanmax_nominal/(granularity * hlat);
    max_memory_multipler = u*v*max_number_elems / (10e3*10e3);
    split_factor = ceil(sqrt(max_memory_multipler));
end
    
if split_factor > 1

    [rowStart, rowEnd, colStart, colEnd] = submatrixLimits(height, split_factor);
    horizon = zeros(v, u, 'like', height);
    col_margin = ceil(spanmax_nominal/hlon);
    row_margin = ceil(spanmax_nominal/hlat);

    tileCount = split_factor^2;
    
    for tileIdx = 1:tileCount

        rowIdxs = max(1, rowStart(tileIdx) - row_margin):min(v, rowEnd(tileIdx) + row_margin);
        colIdxs = max(1, colStart(tileIdx) - col_margin):min(u, colEnd(tileIdx) + col_margin);

        tile_longrid = longrid(rowIdxs, colIdxs);
        tile_latgrid = latgrid(rowIdxs, colIdxs);
        Finterp_height = height(rowIdxs, colIdxs);

        if flag_parallel
            tile_horizon = height2horizon_par(tile_longrid, tile_latgrid, Finterp_height, struct('granularity', granularity,'spanmax', spanmax));
        else
            tile_horizon = height2horizon_vec(tile_longrid, tile_latgrid, Finterp_height, struct('granularity', granularity,'spanmax', spanmax));
        end

        if rowStart(tileIdx) == 1; subrowStart = 1; else subrowStart = 1 + row_margin; end
        if rowEnd(tileIdx) == v; subrowEnd = length(rowIdxs); else subrowEnd = length(rowIdxs) - row_margin; end
        if colStart(tileIdx) == 1; subcolStart = 1; else subcolStart = 1 + col_margin; end
        if colEnd(tileIdx) == u; subcolEnd = length(colIdxs); else subcolEnd = length(colIdxs) - col_margin; end
        subrowIdxs = subrowStart:subrowEnd;
        subcolIdxs = subcolStart:subcolEnd;
        horizon(rowIdxs(subrowIdxs), colIdxs(subcolIdxs)) = tile_horizon(subrowIdxs, subcolIdxs);

        if flag_plot
            if tileIdx == 1
                figure(), 
                hold on, axis equal, 
                set(gca(), 'YDir','normal')
                xlim([1, u])
                xlabel('[px]')
                ylim([1, v])
                ylabel('[px]')
                cb = colorbar();
                cb.Limits = [-90, 90];
            end
            imagesc(colIdxs, rowIdxs, rad2deg(horizon(rowIdxs, colIdxs)))
            drawnow
        end

        fprintf('\n%.1f%%. Tile %d of %d: rows %d - %d / columns %d - %d', 100*tileIdx/tileCount, tileIdx, tileCount, rowIdxs(1), rowIdxs(end), colIdxs(1), colIdxs(end));
    end

else        
    if flag_parallel
        [horizon] = height2horizon_par(longrid, latgrid, height, struct('granularity',granularity,'spanmax', spanmax));
    else
        [horizon] = height2horizon_vec(longrid, latgrid, height, struct('granularity',granularity,'spanmax', spanmax));
    end
end

scaling = pi/2; % this makes the data bound between -1 and 1

try
    tiffCorrect(filepath_dem, filepath_horizon, single(horizon/scaling), nbit_map);
catch
    warning("Couldn't save grayscale image")
end

if flag_plot
    figure(), 
    hold on, axis equal, 
    imagesc(rad2deg(lon), rad2deg(lat), height - mean(Rbody))
    cb = colorbar;
    cb.Label.String = 'Displacement [m]';
    colormap('turbo')
    xlim([rad2deg(dem_limits(1,:))])
    ylim([rad2deg(dem_limits(2,:))])
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    title('Displacement Map')

    figure(), 
    hold on, axis equal, 
    imagesc(rad2deg(lon), rad2deg(lat), rad2deg(horizon))
    xlim([rad2deg(dem_limits(1,:))])
    ylim([rad2deg(dem_limits(2,:))])
    cb = colorbar;
    cb.Label.String = 'Elevation [deg]';
    colormap('turbo')
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    title('Horizon Map')
end

end