function [valTot, valVec, valROI, nThresholded] = extractImageROI(img, centerROI_px, sizeROI_px, offsetROI_px, valThreshold, flag_plot)

if ~exist('flag_plot','var')
    flag_plot = false;
end

if ~iscell(img)
    img = {img};
end

for ix = 1:length(img)

    imgTemp = img{ix};
    if size(centerROI_px, 2) == 1
        bodyPxCenter_temp = centerROI_px; 
    else 
        bodyPxCenter_temp = centerROI_px(:, ix); 
    end
    if length(sizeROI_px) == 1
        bodySize_px_temp = sizeROI_px; 
    else 
        bodySize_px_temp = sizeROI_px(ix); 
    end
    if length(offsetROI_px) == 1
        roiWindowOffset_temp = offsetROI_px; 
    else 
        roiWindowOffset_temp = offsetROI_px(ix); 
    end
    if length(valThreshold) == 1
        dnThreshold_temp = valThreshold; 
    else 
        dnThreshold_temp = valThreshold(ix); 
    end
    roiPxRadius = ceil(bodySize_px_temp/2) + round(roiWindowOffset_temp);
    roiPxVec = [-roiPxRadius:1:roiPxRadius];
    ixsCol = ceil(bodyPxCenter_temp(1)) + roiPxVec;
    ixsRow = ceil(bodyPxCenter_temp(2)) + roiPxVec;
    imgROI = imgTemp(ixsRow, ixsCol);
    ixsThresholded = imgROI>dnThreshold_temp;
    imgROI(~ixsThresholded) = 0;
    nThresholded(ix) = sum(ixsThresholded, 'all');

    valROI{ix} = imgROI;
    valVec(1:nThresholded(ix), ix) = imgROI(ixsThresholded);
    valTot(ix) = sum(valVec(1:nThresholded(ix), ix));

    if flag_plot
        figure('name','roi','units','pixels','Position',[100 100 300 size(imgTemp,2)/size(imgTemp,1)*300])
        set(gca(), 'Position',[0 0 1 1])
        if max(imgROI(:))>255
            imagesc(imgTemp(ixsRow, ixsCol))
            colormap('jet')        
            ylim([0.5, size(imgROI, 1)])
            xlim([0.5, size(imgROI, 2)])
            col = 'white';
        else
            imshow(imgTemp(ixsRow, ixsCol),'InitialMagnification', 'fit');
            col = 'red';
        end
        hold on, axis equal
        [B,~] = bwboundaries(imgROI, 'noholes'); 
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), col, 'LineWidth', 2);
        end                
        scatter(size(imgROI, 1)/2 + 0.5, size(imgROI, 2)/2 + 0.5, 100, '+', 'LineWidth', 2,'MarkerEdgeColor',col);
    end

end

end