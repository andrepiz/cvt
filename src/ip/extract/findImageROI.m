function [uvROI_px, sizeROI_px, dnThreshold, dnThresholdBackground] = ...
    findImageROI(img, fltSize, bodySizeApprox_px, noiseThreshold, method, ...
                phaseAngleApprox, flag_plot)

if ~exist('fltSize','var')
    fltSize = [];
end
if ~exist('bodySizeApprox_px','var')
    bodySizeApprox_px = [];
end
if ~exist('noiseThreshold','var')
    noiseThreshold = 'dynamic';
end
if ~exist('method','var')
    method = 'moments';
end
if ~exist('phaseAngleApprox','var')
    phaseAngleApprox = 0;
end
if ~exist('flag_plot','var')
    flag_plot = false;
end

if ~iscell(img)
    img = {img};
end

phaseAngleScaling = (1 + cos(phaseAngleApprox))/2;

for ix = 1:length(img)

    if isempty(bodySizeApprox_px)
        % Unknown body size, use all the pixels in the image
        nfrac = 1;
    else
        if length(bodySizeApprox_px) == 1
            bodySizeApprox_px_temp = bodySizeApprox_px; 
        else 
            bodySizeApprox_px_temp = bodySizeApprox_px(ix); 
        end    
        % Expected percentage of pixels illuminated by the target
        nfrac = phaseAngleScaling*pi/4*(bodySizeApprox_px_temp).^2/numel(img{ix});
    end

    if isempty(fltSize)
        % Unknown filter size, do not apply a box filter
        fltSize_temp = 1;
    else
        if length(fltSize) == 1
            fltSize_temp = fltSize; 
        else 
            fltSize_temp = fltSize(ix); 
        end    
    end

    % Centroiding
    [uvROIAll_px, dnThreshold(ix), dnThresholdBackground(ix), maskNoiseTemp] = centroidBrightestPixels(img{ix}, fltSize_temp, noiseThreshold, nfrac, method, flag_plot);

    if isempty(uvROIAll_px)
        warning('NO ROI FOUND')
        uvROI_px(:, ix) = [nan; nan];
        sizeROI_px(1, ix) = nan;
        continue
    else
        uvROI_px(:, ix) = uvROIAll_px(:, 1); % taking the brightest ROI
        sizeROI_px(1, ix) = 2*sqrt(1/phaseAngleScaling*numel(img{ix}(~maskNoiseTemp))/pi);
    end

    if flag_plot
        th = linspace(0, 2*pi, 100);
        figure('Name','centroiding_roi'), 
        subplot(1,2,1)
        imagesc(img{ix});
        hold on, axis equal
        scatter(uvROI_px(1, ix), uvROI_px(2, ix), 50, 'r+')
        plot(uvROI_px(1, ix) + sizeROI_px(1, ix)/2 * cos(th), uvROI_px(2, ix) + sizeROI_px(1, ix)/2 * sin(th), 'c--', 'LineWidth', 2);
        ylim([0, size(img{ix}, 1)])
        xlim([0, size(img{ix}, 2)])
        subplot(1,2,2)
        imagesc(img{ix});
        hold on, axis equal
        scatter(uvROI_px(1, ix), uvROI_px(2, ix), 50, 'r+')
        plot(uvROI_px(1, ix) + sizeROI_px(1, ix)/2 * cos(th), uvROI_px(2, ix) + sizeROI_px(1, ix)/2 * sin(th), 'c--', 'LineWidth', 2);
        xlim(uvROI_px(1, ix) + [-1 1]*sizeROI_px(1, ix))
        ylim(uvROI_px(2, ix) + [-1 1]*sizeROI_px(1, ix))
    end
end

end