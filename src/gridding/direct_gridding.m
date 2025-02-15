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
dValsScaled = dValues*dGranularity;

switch chMethod
    case 'sum'

        valsPixelScaledFine = histsum_2d(dCoords, dValsScaled, dLimits, dGranularity);

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
        if bParallelization && length(dValsScaled) > 1e6    % Over about 1M points the parallelization version is not faster
             valsPixelScaledFine = parhistweight_2d(dCoords, dValsScaled, dLimits, dGranularity, ...
                                                    i32Algorithm, dWindow, 1/2, ...
                                                    false, ui8Workers);
        else
             valsPixelScaledFine = histweight_2d(dCoords, dValsScaled, dLimits, dGranularity, ...
                                                 i32Algorithm, dWindow, 1/2, ...
                                                 false);
        end

    case 'interpolation'
        valsPixelScaledFine = quantization(dCoords, dValsScaled, dLimits, dGranularity, 'method', chScheme);

    case 'shiftedsum'
        valsPixelScaledFine = shiftedquantization(dCoords, dValsScaled, dLimits, dGranularity, 'method', 'sum', 'shift', dShift);

    case 'weightedshiftedsum'
        switch chAlgorithm
            case 'gaussian'
                kern_gridding = gaussianKernel(dShift, sqrt(dShift), false);
            otherwise
                error('Gridding filter not supported')
        end
        valsPixelScaledFine = shiftedquantization(dCoords, dValsScaled, dLimits, dGranularity, 'method', 'sum', 'shift', dShift, 'weight', kern_gridding);

    otherwise
        error('Gridding method not recognized')
end

if dGranularity == 1
    valsPixelScaled = valsPixelScaledFine;
else
    if bAntialiasing
       kern_antialiasing = gaussianKernel(dGranularity, sqrt(dGranularity), false);
    else
       kern_antialiasing = [];
    end
    valsPixelScaled = downsamplingreconstruction(valsPixelScaledFine, dGranularity, kern_antialiasing, chFilter, dGranularity);
    valsPixelScaled(valsPixelScaled<0) = 0;
end

% removing nans
ixsNan = isnan(valsPixelScaled);
valsPixelScaled(ixsNan) = 0;

% Reapply scaling factor
mat = valsPixelScaled/dGranularity;
maskValid = ~ixsNan;

end