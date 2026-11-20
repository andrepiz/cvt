function [uvROI_px, sideROI_px, dnThreshold, dnThresholdBackground, nPxLit, nPxDisk] = ...
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

nimg = length(img);

phaseAngleScaling = (1 + cos(phaseAngleApprox))/2;
if length(phaseAngleScaling) == 1
    phaseAngleScaling = repmat(phaseAngleScaling, 1, nimg);
end

flag_disp = isempty(gcp('nocreate')) && nimg > 1;
if ~flag_disp && nimg > 1; fprintf('\nCentroiding %d images...', nimg); end

parfor ix = 1:nimg

    if nimg > 1 && flag_disp
        fprintf('\nCentroiding image %d out of %d', ix, nimg);
    end

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
        nfrac = phaseAngleScaling(ix)*pi/4*(bodySizeApprox_px_temp).^2/numel(img{ix});
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
    [uvROIAll_px, dnThreshold(ix), dnThresholdBackground(ix), maskNoiseTemp, valLocTemp] = centroidBrightestPixels(img{ix}, fltSize_temp, noiseThreshold, nfrac, method, false);
    %noiseThreshold = dnThreshold(ix); % using previous noise threshold (images provided as sequence)

    if isempty(uvROIAll_px)
        warning('NO ROI FOUND')
        uvROI_px(:, ix) = [nan; nan];
        nPxLit(1, ix) = nan;
        nPxDisk(1, ix) = nan;
        sideROI_px(1, ix) = nan;
        continue
    else
        uvROI_px(:, ix) = uvROIAll_px(:, 1); % taking the brightest ROI
        nPxLit(1, ix) = numel(img{ix}(~maskNoiseTemp));
        nPxDisk(1, ix) = nPxLit(1, ix)./phaseAngleScaling(ix);
        sideROI_px(1, ix) = 2*sqrt(nPxDisk(1, ix)/pi);
    end

end

if flag_plot

    for ix = 1:nimg
        th = linspace(0, 2*pi, 100);
        figure('Name','centroiding_roi'), 
        subplot(1,2,1)
        imagesc(img{ix});
        hold on, axis equal
        scatter(uvROI_px(1, ix), uvROI_px(2, ix), 50, 'r+')
        plot(uvROI_px(1, ix) + sideROI_px(1, ix)/2 * cos(th), uvROI_px(2, ix) + sideROI_px(1, ix)/2 * sin(th), 'c--', 'LineWidth', 2);
        ylim([0, size(img{ix}, 1)])
        xlim([0, size(img{ix}, 2)])
        subplot(1,2,2)
        imagesc(img{ix});
        hold on, axis equal
        scatter(uvROI_px(1, ix), uvROI_px(2, ix), 50, 'r+')
        plot(uvROI_px(1, ix) + sideROI_px(1, ix)/2 * cos(th), uvROI_px(2, ix) + sideROI_px(1, ix)/2 * sin(th), 'c--', 'LineWidth', 2);
        xlim(uvROI_px(1, ix) + [-1 1]*sideROI_px(1, ix))
        ylim(uvROI_px(2, ix) + [-1 1]*sideROI_px(1, ix))
    end
end

end