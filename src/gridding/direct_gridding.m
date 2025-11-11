function [dMat, dMaskValid] = direct_gridding(dCoordsRC, dVals, dLimsRC,...
                                   bParallelization, ui8Workers, ...
                                  chMethod, dWindow, chAlgorithm, chScheme, dShift, ...
                                  dGranularity, bAntialiasing, chFilter, dSigma)
% DIRECT_GRIDDING Wrapper function to perform direct gridding rasterization using 
%different methods.
%
% INPUTS:
%   dCoordsRC     - 2xN array of decimal rows and columns coordinates
%   dValues       - 1xN vector of values
%   dLimsRC       - 2x2 array of [row_min row_max; col_min col_max] limits
%                   (default: tight bounds on dCoordsRC)
%
% OUTPUTS:
%   dMat          - 2D grid of gridded values
%   dMaskValid    - Logical mask of valid (non-NaN) pixels

arguments
    dCoordsRC           (:, :) double {ismatrix}
    dVals               (1, :) double {isvector}
    dLimsRC             (:, 2) double {ismatrix} = [floor(min(dCoordsRC,[],2)), 1 + ceil(max(dCoordsRC,[],2))];
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
    dSigma              (1, 1) double {isscalar}  = 0.5
end

% Scaling of values
if dGranularity == 1 || strcmp(chMethod, 'aggregation')
    dValuesScaled = dVals;
else
    dValuesScaled = dVals*dGranularity;
end

switch chMethod
    case 'aggregation'
        valsPixelScaledFine = aggregation_2d(dCoordsRC, dValuesScaled, dLimsRC, chAlgorithm, dGranularity);

    case 'sum'
        valsPixelScaledFine = histsum_2d(dCoordsRC, dValuesScaled, dLimsRC, dGranularity);

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
            otherwise
                error('Weighted sum algorithm not recognized')
        end
        if bParallelization && length(dValuesScaled) > 1e6    % Over about 1M points the parallelization version is not faster
             valsPixelScaledFine = parhistweight_2d(dCoordsRC, dValuesScaled, dLimsRC, dGranularity, ...
                                                    i32Algorithm, dWindow, dSigma, ...
                                                    false, ui8Workers);
        else
             valsPixelScaledFine = histweight_2d(dCoordsRC, dValuesScaled, dLimsRC, dGranularity, ...
                                                 i32Algorithm, dWindow, dSigma, ...
                                                 false);
        end

    case 'interpolation'
        valsPixelScaledFine = quantization(dCoordsRC, dValuesScaled, dLimsRC, dGranularity, 'method', chScheme);

    case 'shiftedsum'
        valsPixelScaledFine = shiftedquantization(dCoordsRC, dValuesScaled, dLimsRC, dGranularity, 'method', 'sum', 'shift', dShift);

    case 'weightedshiftedsum'
        switch chAlgorithm
            case 'gaussian'
                kern_gridding = gaussianKernel(dShift, sqrt(dShift), false);
            otherwise
                error('Gridding filter not supported')
        end
        valsPixelScaledFine = shiftedquantization(dCoordsRC, dValuesScaled, dLimsRC, dGranularity, 'method', 'sum', 'shift', dShift, 'weight', kern_gridding);

    otherwise
        error('Gridding method not recognized')
end

if dGranularity == 1
    dMat = valsPixelScaledFine;
else
    % Create reconstruction anti-aliasing kernel if needed
    if bAntialiasing
       kern_antialiasing = gaussianKernel(dGranularity, sqrt(dGranularity), false);
    else
       kern_antialiasing = [];
    end
    % Reconstruct the matrix to its original size
    dMat = downsamplingreconstruction(valsPixelScaledFine, dGranularity, kern_antialiasing, chFilter, dGranularity);
    dMat(dMat<0) = 0;
end

if ~strcmp(chMethod, 'aggregation')
    % Reapply scaling factor
    dMat = dMat/dGranularity;
end

% setting nans to 0
ixsNan = isnan(dMat);
dMaskValid = ~ixsNan;
dMat(ixsNan) = 0;

end