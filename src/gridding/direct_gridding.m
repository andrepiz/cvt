function [mat, maskValid] = direct_gridding(dCoords, dValues, dLimits,...
                                   bParallelization, ui8Workers, ...
                                  chMethod, dWindow, chAlgorithm, chScheme, dShift, ...
                                  dGranularity, bAntialiasing, chFilter)

arguments
    dCoords             (:, :) double {ismatrix}
    dValues             (1, :) double {isvector}
    dLimits             (:, 2) double {ismatrix} = [floor(min(dCoords,[],2)), 1 + ceil(max(dCoords,[],2))];
    bParallelization    (1, 1) logical           = false
    ui8Workers          (1, 1) uint8 {isscalar}  = 4   
    chMethod                   char              = 'weightedsum'
    dWindow             (1, 1) double {isscalar}  = 1
    chAlgorithm                char              = 'gaussian'   
    chScheme                   char              = 'linear'
    dShift              (1, 1) double {isscalar}  = 1
    dGranularity        (1, 1) double {isscalar}  = 1
    bAntialiasing       (1, 1) logical           = false
    chFilter                   char              = 'gaussian'
end


% Scaling of values
vals_scaled = dValues*dGranularity;

switch chMethod
    case 'sum'
        [valsPixel_scaled_fine, counts, edges] = quantization(dCoords, vals_scaled, dLimits, dGranularity, 'method', 'sum');

    case 'weightedsum'
        switch chAlgorithm
            case 'invsquared'
                i32Algorithm = 0;
            case 'diff'
                i32Algorithm = 1;
            case 'area'
                i32Algorithm = 2;
            case 'gaussian'
                i32Algorithm = 3;
        end
        if bParallelization
            [valsPixel_scaled_fine, counts, edges] = parhistweight_2d(dCoords, vals_scaled, dLimits, dGranularity, i32Algorithm, ui8Workers, 1/3, dWindow);
        else
            [valsPixel_scaled_fine, counts] = histweight_2d(dCoords, vals_scaled, dLimits, dGranularity, i32Algorithm, false, false, false, 1/3, dWindow);
        end

    case 'interpolation'
        [valsPixel_scaled_fine, counts, edges] = quantization(dCoords, vals_scaled, dLimits, dGranularity, 'method', chScheme);

    case 'shiftedsum'
        [valsPixel_scaled_fine, counts, edges] = shiftedquantization(dCoords, vals_scaled, dLimits, dGranularity, 'method', 'sum', 'shift', dShift);

    case 'weightedshiftedsum'
        switch chAlgorithm
            case 'gaussian'
                kern_gridding = gaussianKernel(dShift, sqrt(dShift), false);
            otherwise
                error('Gridding filter not supported')
        end
        [valsPixel_scaled_fine, counts, edges] = shiftedquantization(dCoords, vals_scaled, dLimits, dGranularity, 'method', 'sum', 'shift', dShift, 'weight', kern_gridding);

    otherwise
        error('Gridding method not recognized')
end

if dGranularity == 1
    valsPixel_scaled = valsPixel_scaled_fine;
else
    if bAntialiasing
       kern_antialiasing = gaussianKernel(dGranularity, sqrt(dGranularity), false);
    else
       kern_antialiasing = [];
    end
    valsPixel_scaled = downsamplingreconstruction(valsPixel_scaled_fine, dGranularity, kern_antialiasing, chFilter, dGranularity);
    valsPixel_scaled(valsPixel_scaled<0) = 0;
end

% removing nans
ixsNan = isnan(valsPixel_scaled);
valsPixel_scaled(ixsNan) = 0;

% Reapply scaling factor
mat = valsPixel_scaled/dGranularity;
maskValid = ~ixsNan;

end