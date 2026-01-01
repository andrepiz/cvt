function fh = addZoomInsetToFigure(fh, posCrop, posInset, lw)

[xCrop, yCrop, widthCrop, heightCrop] = deal(posCrop(1), posCrop(2), posCrop(3), posCrop(4));
[xInset, yInset, widthInset, heightInset] = deal(posInset(1), posInset(2), posInset(3), posInset(4));
axMain = fh.CurrentAxes;
dataInset = allchild(axMain);

line(axMain, [xCrop, xInset], [yCrop, yInset], 'Color', 'r', 'LineWidth', lw);
line(axMain, [xCrop + widthCrop, xInset + widthInset], [yCrop, yInset], 'Color', 'r', 'LineWidth', lw);
line(axMain, [xCrop + widthCrop, xInset + widthInset], [yCrop + heightCrop, yInset + heightInset], 'Color', 'r', 'LineWidth', lw);
line(axMain, [xCrop, xInset], [yCrop + heightCrop, yInset + heightInset], 'Color', 'r', 'LineWidth', lw);
rectangle(axMain, 'Position', [xCrop, yCrop, widthCrop, heightCrop], 'EdgeColor', 'r', 'LineWidth', lw);
rectangle(axMain, 'Position', [xInset, yInset, widthInset, heightInset], 'EdgeColor', 'r', 'LineWidth', 2*lw);

insetPos = insetAxesDataToFigurePixels(axMain, [xInset, yInset, widthInset, heightInset]);

[xInsetPos, yInsetPos, widthInsetPos, heightInsetPos] = deal(insetPos(1), insetPos(2), insetPos(3), insetPos(4));
axInset = axes('units','pixel','Position', [xInsetPos, yInsetPos, widthInsetPos, heightInsetPos]);
copyobj(dataInset, axInset);
hold on, axis equal, axis off
set(axInset,'XLim',[axMain.XLim],'YLim',[axMain.YLim],'YDir','reverse')
xlim(axInset, xCrop + [0, widthCrop])
ylim(axInset, yCrop + [0, heightCrop])

end