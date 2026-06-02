function I_lin = srgb2lin(I_n)
% SRGB2LIN  Convert sRGB-encoded image to linear light values.
%
%   I_lin = srgb2lin(I_n) applies the inverse sRGB transfer function
%   (IEC 61966-2-1) to a normalized input image I_n in [0, 1],
%   returning linear light values I_lin in [0, 1].
%
%   Input:
%     I_n   - sRGB-encoded image, double, values in [0, 1]
%   Output:
%     I_lin - linearized image (scene-linear light), same size as I_n
%
%   Note: negative values or values > 1 will not be clamped.
%         Use max(0, min(1, I_n)) beforehand if needed.

if min(I_n(:)) < 0 || max(I_n(:)) > 1
    error('Input image should be normalized between 0 and 1.')
end

a = 0.055;
thr = 0.04045;
I_lin = zeros(size(I_n));

mask = I_n <= thr;
I_lin(mask)  = I_n(mask) / 12.92;
I_lin(~mask) = ((I_n(~mask) + a) / (1 + a)).^2.4;

end