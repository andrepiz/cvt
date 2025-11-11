function ixsUnoccluded = depthbuffering(dCoordsRC, dDepth, dLimitsRC, dCoordsBiasedRC, dDepthBiased, dDepthDiffThreshold, dGranularity)
% Depth buffer is created using original coords and depths, but depth check
% is done using modified coords and depths. 
% In case bias is along normals: points at nadir will have same coords but
% smaller depth than the one in the buffer, so they will be visible; points
% off-nadir will have coords shifted and possibly larger depths

% CREATE DEPTH BUFFER WITH ORIGINAL COORDS-DEPTHS
coords = dCoordsRC*dGranularity;
vals = dDepth;
lims  = [(dLimitsRC(:,1) - 1)*dGranularity, dLimitsRC(:, 2)*dGranularity];

% pts_bounded
coords_shifted = coords - lims(:, 1);
row_idx = ceil(coords_shifted(1, :));  % Point at [0.3, 0.7] will snap to row 1 and col 1
col_idx = ceil(coords_shifted(2, :));

% Restrict to submatrix
nrows = max(row_idx);
ncols = max(col_idx);
dBins = zeros(nrows, ncols);
nbins = numel(dBins);

% Linear indexing and aggregation into depth buffer
idx = sub2ind([nrows, ncols], row_idx, col_idx);
binsVec = accumarray(idx', vals', [nbins, 1], @min, NaN);
depthBuffer = reshape(binsVec, [nrows, ncols]);

if ~isempty(dCoordsBiasedRC)
    coordsBiased = dCoordsBiasedRC*dGranularity;
    coordsBiased_shifted = coordsBiased - lims(:, 1);
    rowBiased_idx = ceil(coordsBiased_shifted(1, :));  % Point at [0.3, 0.7] will snap to row 1 and col 1
    colBiased_idx = ceil(coordsBiased_shifted(2, :));
    
    ixsInside = rowBiased_idx >= 1 & rowBiased_idx <= nrows & ...
                colBiased_idx >= 1 & colBiased_idx <= ncols;
    
    % Linear indexing
    idxBiased = sub2ind([nrows, ncols], rowBiased_idx(ixsInside), colBiased_idx(ixsInside));

    % Smooth original depthBuffer as we need to read from different rows
    % and cols
    depthBuffer = fillmissing2(depthBuffer, "movmean",3);
    % figure(), imagesc(depthBuffer - reshape(binsVec, [nrows, ncols]))
else
    idxBiased = idx;
    ixsInside = true(size(dDepth));
    dDepthBiased = dDepth;
end

depthBiasedBufferCheckInside = depthBuffer(idxBiased);

depthDiffInside = dDepthBiased(ixsInside) - depthBiasedBufferCheckInside;

ixsUnoccluded = true(size(dDepthBiased));
if numel(dDepthDiffThreshold) == 1
    dDepthDiffThreshold = dDepthDiffThreshold*ones(size(vals));
end
ixsUnoccluded(ixsInside) = depthDiffInside <= dDepthDiffThreshold(ixsInside);
%ixsUnoccluded(isnan(depthDiffBufferCheck)) = true;    
    
flag_debug = false;
if flag_debug

    ixsPlot = round(linspace(1,length(vals), 1e5));
    coordsInside = coords(:, ixsInside);
    ixsPlotInside = round(linspace(1,size(coordsInside, 2), 1e5));

    fh1 = figure();
           
    ax0 = subplot(2,2,1);
    grid on, hold on, axis equal
    scatter(gca(), dCoordsRC(2,ixsPlot), dCoordsRC(1,ixsPlot), [], vals(ixsPlot))
    ax0.YDir = 'reverse';
    colorbar()
    xlabel('u [px]')
    ylabel('v [px]')   

    ax1 = subplot(2,2,2);
    grid on, hold on, axis equal
    imagesc(depthBuffer)
    ax1.YDir = 'reverse';
    colorbar()
    xlabel('u [px]')
    ylabel('v [px]')            

    ax2 = subplot(2,2,3);
    grid on, hold on, axis equal
    grid minor
    scatter(gca(), coordsInside(2,ixsPlotInside), coordsInside(1,ixsPlotInside), [], depthDiffInside(ixsPlotInside))
    axis equal
    ax2.YDir = 'reverse';
    col = colorbar;
    col.Label.String = '[-]';
    xlabel('u [px]')
    ylabel('v [px]')            
    xlim(lims(2,:))
    ylim(lims(1,:))

    try
    ax3 = subplot(2,2,4);
    grid on, hold on, axis equal
    grid minor
    scatter(gca(), coords(2,ixsPlot), coords(1,ixsPlot), [], dDepthDiffThreshold(ixsPlot))
    axis equal
    ax3.YDir = 'reverse';
    col = colorbar;
    col.Label.String = '[-]';
    xlabel('u [px]')
    ylabel('v [px]')            
    xlim(lims(2,:))
    ylim(lims(1,:))
    catch 
    end
    
    figure()
    grid on, hold on, axis equal
    grid minor
    scatter(gca(), dCoordsRC(2,~ixsUnoccluded), dCoordsRC(1,~ixsUnoccluded), 'red')
    axis equal
    set(gca(),'YDir', 'reverse');
    col = colorbar;
    col.Label.String = '[-]';
    xlabel('u [px]')
    ylabel('v [px]')            
    xlim(dLimitsRC(2,:))
    ylim(dLimitsRC(1,:))
end

end