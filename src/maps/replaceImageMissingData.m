function I_out = replaceImageMissingData(I, threshold, resize_factor, kernel_size)

I_high = double(I);
sca = max(I_high(:));

I_low = imresize(I_high/sca, 1/resize_factor, {@(x) validKernel(x, threshold/sca), kernel_size});
try
    I_out = substituteLowResNoResizeFast(I_high, I_low*sca, threshold);
catch
    I_out = substituteLowResNoResizeSafe(I_high, I_low*sca, threshold);
end

end

function out = validKernel(x, threshold)

out = zeros(size(x));
out(x >= threshold) = 1;
out = out/sum(out(:));

end

function high_res_filled = substituteLowResNoResizeFast(high_res, low_res, threshold)

    [H_high, W_high] = size(high_res);
    [H_low, W_low] = size(low_res);
    
    mask = high_res < threshold;
    [rows, cols] = find(mask);
    
    % Vectorized coordinate mapping
    i_low = max(1, min(H_low, floor((rows-1) * H_low / H_high) + 1));
    j_low = max(1, min(W_low, floor((cols-1) * W_low / W_high) + 1));
    
    high_res_filled = high_res;
    high_res_filled(sub2ind(size(high_res), rows, cols)) = low_res(sub2ind(size(low_res), i_low, j_low));
    
    num_replaced = length(rows);
    fprintf('Vectorized: Replaced %d pixels (%.2f%%)\n', num_replaced, 100*num_replaced/numel(high_res));
end

function high_res_filled = substituteLowResNoResizeSafe(high_res, low_res, threshold)
    [H_high, W_high] = size(high_res);
    [H_low, W_low] = size(low_res);
    
    high_res_filled = high_res;
    num_replaced = 0;
    
    fprintf('Processing row-by-row...\n');
    
    % Precompute scaling factors
    scale_row = H_low / H_high;
    scale_col = W_low / W_high;
    
    tic;
    for ii = 1:H_high
        % Process one row at a time
        row_mask = high_res(ii, :) < threshold;
        if any(row_mask)
            cols_to_fill = find(row_mask);
            num_replaced = num_replaced + length(cols_to_fill);
            
            % Map coordinates for this row only
            i_low = max(1, min(H_low, floor((ii-1) * scale_row) + 1));
            j_low = max(1, min(W_low, floor((cols_to_fill-1) * scale_col) + 1));
            
            % Fill this row segment
            high_res_filled(ii, cols_to_fill) = low_res(i_low, j_low);
        end
        
        if mod(ii, 100) == 0
            fprintf('Row %d/%d (%.1f%%)\n', ii, H_high, 100*ii/H_high);
        end
    end
    
    elapsed = toc;
    fprintf('Safe: Replaced %d pixels (%.2f%%) in %.2f s\n', ...
            num_replaced, 100*num_replaced/(H_high*W_high), elapsed);
end

