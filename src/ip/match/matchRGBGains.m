function [matched, gains, info] = matchRGBGains( ...
    image, template, maskImage, maskTemplate)
%MATCHRGBGAINS Match RGB distributions without pixel correspondence.
%
% Inputs:
%   image, template         sRGB images, or image file paths.
%   maskImage, maskTemplate Optional logical masks, one per image.
%                           Images can have different dimensions.
%
% Outputs:
%   matched  Corrected sRGB image, double in [0,1].
%   gains    [gainR gainG gainB], applied in linear RGB.
%   info     Optimization diagnostics.
%
% Assumes both inputs use the sRGB color space.
% Requires Image Processing Toolbox.

    if ischar(image) || isstring(image)
        image = imread(image);
    end
    if ischar(template) || isstring(template)
        template = imread(template);
    end

    image = im2double(image);
    template = im2double(template);

    assert(size(image,3) == 3 && size(template,3) == 3, ...
        'Inputs must be RGB images.');

    if nargin < 3 || isempty(maskImage)
        maskImage = true(size(image,1), size(image,2));
    end
    if nargin < 4 || isempty(maskTemplate)
        maskTemplate = true(size(template,1), size(template,2));
    end

    assert(isequal(size(maskImage), size(image,[1 2])), ...
        'maskImage must match the image dimensions.');
    assert(isequal(size(maskTemplate), size(template,[1 2])), ...
        'maskTemplate must match the template dimensions.');

    % Remove very dark pixels and pixels clipped in any channel.
    maskImage = logical(maskImage) ...
        & all(isfinite(image),3) ...
        & max(image,[],3) > 0.03 ...
        & max(image,[],3) < 0.98;

    maskTemplate = logical(maskTemplate) ...
        & all(isfinite(template),3) ...
        & max(template,[],3) > 0.03 ...
        & max(template,[],3) < 0.98;

    assert(nnz(maskImage) >= 100 && nnz(maskTemplate) >= 100, ...
        'Each mask must contain at least 100 usable pixels.');

    sourceLinear = rgb2lin(image);
    targetLinear = rgb2lin(template);

    % Match distributions rather than spatially corresponding samples.
    percentiles = 10:10:90;
    sourceQ = zeros(numel(percentiles),3);
    targetQ = zeros(numel(percentiles),3);

    for channel = 1:3
        a = sourceLinear(:,:,channel);
        b = targetLinear(:,:,channel);

        sourceQ(:,channel) = reshape( ...
            prctile(a(maskImage),percentiles), [], 1);
        targetQ(:,channel) = reshape( ...
            prctile(b(maskTemplate),percentiles), [], 1);
    end

    % Log gains ensure positivity. Log residuals measure relative error.
    floorValue = 1e-6;
    logSource = log(max(sourceQ,floorValue));
    logTarget = log(max(targetQ,floorValue));

    options = optimset( ...
        'Display','off', ...
        'MaxIter',1000, ...
        'MaxFunEvals',3000, ...
        'TolX',1e-8, ...
        'TolFun',1e-10);

    [logGains, loss, exitflag, output] = ...
        fminsearch(@objective, zeros(1,3), options);

    gains = exp(logGains);

    correctedLinear = sourceLinear .* reshape(gains,1,1,3);

    % Measure clipping before limiting the output for display.
    clipped = any(correctedLinear > 1,3);
    info.clippedFraction = nnz(clipped & maskImage) / nnz(maskImage);

    matched = lin2rgb(min(max(correctedLinear,0),1));

    info.initialLoss = objective(zeros(1,3));
    info.finalLoss = loss;
    info.exitflag = exitflag;
    info.iterations = output.iterations;
    info.percentiles = percentiles;
    info.sourceMask = maskImage;
    info.templateMask = maskTemplate;

    function value = objective(x)
        residual = logSource + reshape(x,1,3) - logTarget;

        % Pseudo-Huber loss reduces sensitivity to mismatched percentiles.
        delta = 0.15;
        penalty = delta^2 * ...
            (sqrt(1 + (residual/delta).^2) - 1);

        value = mean(penalty(:));
    end
end