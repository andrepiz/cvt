function [bins, counts, edges] = quantization(coords, values, limits, granularity, params)
% Group and sum scattered data points into fixed-size quantiles with specified granularity.
% Each point value is kept in the nearest-neighbour sector or
% interpolated with a different params.method

arguments
    coords double
    values double
    limits double = [floor(min(coords,[],2)), 1 + ceil(max(coords,[],2))]
    granularity (1,1) double = 1
    params.method char = 'gaussian'
    params.flag_progress (1,1) logical = false
end

coords_scaled = coords*granularity;
nx = limits(1,2)*granularity;
ny = limits(2,2)*granularity;

% Partitioning
xsubedges = {linspace(limits(1,1) - 1, nx/granularity, nx + 1)};
ysubedges = {linspace(limits(2,1) - 1, ny/granularity, ny + 1)};
c = 1;
while any(cellfun(@(x) length(x), xsubedges) > 1024) || any(cellfun(@(x) length(x), ysubedges) > 1024)
    % Partition into submatrices
    xstep = floor(nx/(c + 1));
    ystep = floor(ny/(c + 1));
    c = c + 1;
    for ix = 1:c
        xsubedges{ix} = linspace(limits(1,1) - 1 + (ix-1)*xstep, ix*xstep, xstep + 1);
        ysubedges{ix} = linspace(limits(2,1) - 1 + (ix-1)*ystep, ix*ystep, ystep + 1);
    end
end
counts = zeros(nx, ny);
bins = zeros(nx, ny);
edges = {unique([xsubedges{:}]), unique([ysubedges{:}])};
c = 0;

% Binning
for ii = 1:length(xsubedges)
    for jj = 1:length(ysubedges)

        c = c + 1;
        if mod(c, round((length(xsubedges)*length(ysubedges))/10)) == 0 && params.flag_progress
            disp(['     ',num2str(round(1e2*(c)/(length(xsubedges)*length(ysubedges)))),'%'])
        end
        
        subixs = coords_scaled(1,:) > min(xsubedges{ii}) & coords_scaled(1,:) <= max(xsubedges{ii}) & ...
                 coords_scaled(2,:) > min(ysubedges{jj}) & coords_scaled(2,:) <= max(ysubedges{jj});
        subcoords_scaled = coords_scaled(:, subixs);
        subvalues = values(:, subixs);
        if isempty(subcoords_scaled)
            continue
        end

        [subcounts, ~, ~, idx_x, idx_y] = histcounts2(subcoords_scaled(1,:), subcoords_scaled(2,:), xsubedges{ii}, ysubedges{jj});
        subbins = zeros(size(subcounts));
        switch params.method
            case 'sum'
                subcounts_check = zeros(size(subcounts));
                for ix = 1:length(subvalues)
                    subcounts_check(idx_x(ix), idx_y(ix)) = subcounts_check(idx_x(ix), idx_y(ix)) + 1;
                    subbins(idx_x(ix), idx_y(ix)) = subbins(idx_x(ix), idx_y(ix)) + subvalues(ix);
                end
                if sum(subcounts_check - subcounts,'all') ~= 0 
                    error('Check failed')
                end

            otherwise 
                xsubedge_temp = xsubedges{ii};
                ysubedge_temp = ysubedges{jj};
                xv = (xsubedge_temp(1:end - 1) + xsubedge_temp(2:end))/2;
                yv = (ysubedge_temp(1:end - 1) + ysubedge_temp(2:end))/2;
                [xvg, yvg] = meshgrid(xv, yv);
                subbins = griddata(coords_scaled(1,:), coords_scaled(2,:), values, xvg, yvg, params.method)';
                if isempty(subbins)
                   subbins = zeros(size(subcounts));
                end
        end

        subbinsixs = {(ii-1)*xstep + 1:ii*xstep, (jj-1)*ystep + 1:jj*ystep};
        counts(subbinsixs{:}) = subcounts;
        bins(subbinsixs{:}) = subbins;
    end

end

if ~strcmp(params.method,'sum')
    bins = fillmissing(bins, params.method, 'MaxGap', 1);
    bins = fillmissing(bins, 'constant', 0);
end

end

