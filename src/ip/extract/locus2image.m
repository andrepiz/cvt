function [img_boundary, img_fill] = locus2image(u, v, res_px)
    % locus2image - creates line and filled binary images from locus of points
    %
    % Inputs:
    %   UV - Nx2 array of (U,V) coordinates defining the polygon vertices
    %   img_res - [height, width] of the output pixel image
    %
    % Outputs:
    %   img_boundary - binary image with polygon edges drawn (1 on edges, 0 elsewhere)
    %   img_shape - binary image with polygon interior filled (1 inside, 0 outside)
    %
    % Note: Coordinates are assumed to be in continuous space;
    %       they are scaled and shifted to fit the image resolution.

    W = res_px(1);
    H = res_px(2);
    
    % Clip input coordinates inside image
    u = min(max(round(u),1), W);
    v = min(max(round(v),1), H);
    
    % Initialize images
    img_boundary = false(H,W);
    
    % Draw edges with Bresenham line algorithm
    N = length(u);
    for k = 1:N
        r1 = v(k);
        c1 = u(k);
        if k < N
            r2 = v(k+1);
            c2 = u(k+1);
        else
            r2 = v(1);
            c2 = u(1);
        end
        
        [rr, cc] = bresenham_line(r1, c1, r2, c2);
        idx = rr >= 1 & rr <= H & cc >= 1 & cc <= W;
        img_boundary(sub2ind([H,W], rr(idx), cc(idx))) = true;
    end
    
    % Create filled polygon image using poly2mask
    img_fill = poly2mask(u, v, H, W);

end

function [rr, cc] = bresenham_line(r1, c1, r2, c2)
    % Bresenham's line algorithm in image coordinates (row, col)
    % Returns the row and column indices of pixels on the line segment

    % Adapted from standard Bresenham algorithm for image coordinates

    dr = abs(r2 - r1);
    dc = abs(c2 - c1);

    if c1 < c2
        step_c = 1;
    else
        step_c = -1;
    end

    if r1 < r2
        step_r = 1;
    else
        step_r = -1;
    end

    rr = [];
    cc = [];

    if dc > dr
        % Iterate over columns
        err = dc / 2;
        r = r1;
        for c = c1:step_c:c2
            rr(end+1) = r; %#ok<AGROW>
            cc(end+1) = c; %#ok<AGROW>
            err = err - dr;
            if err < 0
                r = r + step_r;
                err = err + dc;
            end
        end
    else
        % Iterate over rows
        err = dr / 2;
        c = c1;
        for r = r1:step_r:r2
            rr(end+1) = r; %#ok<AGROW>
            cc(end+1) = c; %#ok<AGROW>
            err = err - dc;
            if err < 0
                c = c + step_c;
                err = err + dr;
            end
        end
    end
end
