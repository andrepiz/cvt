function [valTot, valVec, valROI, nThresholded] = extractImageROI(img, centerROI_px, sideROI_px, offsetROI_px, valThreshold, flag_plot)

if ~exist('flag_plot','var')
    flag_plot = false;
end

if ~iscell(img)
    img = {img};
end

nimg = length(img);

valROI = cell(1, nimg);
valVec = [];
valTot = nan(1, nimg);
nThresholded = nan(1, nimg);
for ix = 1:nimg

    imgTemp = img{ix};
    if size(centerROI_px, 2) == 1
        bodyCenter_px_temp = centerROI_px; 
    else 
        bodyCenter_px_temp = centerROI_px(:, ix); 
    end
    if length(sideROI_px) == 1
        bodySize_px_temp = sideROI_px; 
    else 
        bodySize_px_temp = sideROI_px(ix); 
    end
    if length(offsetROI_px) == 1
        roiOffset_px_temp = offsetROI_px; 
    else 
        roiOffset_px_temp = offsetROI_px(ix); 
    end
    if length(valThreshold) == 1
        dnThreshold_temp = valThreshold; 
    else 
        dnThreshold_temp = valThreshold(ix); 
    end
    roiHalfSize_px_temp = ceil(bodySize_px_temp/2) + round(roiOffset_px_temp);
    ixsColMin = max(1, min(size(imgTemp, 2), ceil(bodyCenter_px_temp(1)) - roiHalfSize_px_temp));
    ixsColMax = max(1, min(size(imgTemp, 2), ceil(bodyCenter_px_temp(1)) + roiHalfSize_px_temp));
    ixsRowMin = max(1, min(size(imgTemp, 1), ceil(bodyCenter_px_temp(2)) - roiHalfSize_px_temp));
    ixsRowMax = max(1, min(size(imgTemp, 1), ceil(bodyCenter_px_temp(2)) + roiHalfSize_px_temp));
    ixsRow = ixsRowMin:ixsRowMax;
    ixsCol = ixsColMin:ixsColMax;
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
        % if max(imgROI(:))>255
            imagesc(imgROI)
            colormap('jet')        
            ylim([1, size(imgROI, 1)])
            xlim([1, size(imgROI, 2)])
            col = 'white';
        % else
        %     imshow(imgTemp(ixsRow, ixsCol),'InitialMagnification', 'fit');
        %     col = 'red';
        % end
        hold on, axis equal
        [B,~] = bwboundaries(imgROI, 'noholes'); 
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), col, 'LineWidth', 2);
        end                
        scatter(size(imgROI, 1)/2+0.5, size(imgROI, 2)/2+0.5, 100, '+', 'LineWidth', 2,'MarkerEdgeColor',col);
    end

end

end