function insetFigurePixels = insetAxesDataToFigurePixels(ax, insetAxesData)

    % Outer axes box in figure pixels
    axPx = getpixelposition(ax, true);   % [axX axY axW axH]

    % Use axis limits as data range
    xl = xlim(ax);
    yl = ylim(ax);

    % Data rectangle
    x1 = insetAxesData(1);
    y1 = insetAxesData(2);
    x2 = x1 + insetAxesData(3);
    y2 = y1 + insetAxesData(4);

    % Map X (same for both YDir)
    px1 = axPx(1) + (x1 - xl(1)) / diff(xl) * axPx(3);
    px2 = axPx(1) + (x2 - xl(1)) / diff(xl) * axPx(3);

    % Map Y with YDir handling
    if strcmp(get(ax,'YDir'),'normal')
        py1 = axPx(2) + (y1 - yl(1)) / diff(yl) * axPx(4);
        py2 = axPx(2) + (y2 - yl(1)) / diff(yl) * axPx(4);
        yInsetPos      = py1;
        heightInsetPos = py2 - py1;
    else  % 'reverse'
        py1 = axPx(2) + (yl(2) - y1) / diff(yl) * axPx(4);
        py2 = axPx(2) + (yl(2) - y2) / diff(yl) * axPx(4);
        yInsetPos      = py2;
        heightInsetPos = py1 - py2;
    end

    xInsetPos     = px1;
    widthInsetPos = px2 - px1;

    insetFigurePixels = [xInsetPos, yInsetPos, widthInsetPos, heightInsetPos];
end
