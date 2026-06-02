function I_enc = lin2srgb(I_lin)
% LIN2SRGB  Apply sRGB transfer function to linear light values.
%
%   I_enc = lin2srgb(I_lin) encodes a linear image I_lin in [0, 1]
%   into sRGB, per IEC 61966-2-1.
%
%   Input:
%     I_lin - linear light image, double, values in [0, 1]
%   Output:
%     I_enc - sRGB-encoded image, same size as I_lin

if min(I_lin(:)) < 0 || max(I_lin(:)) > 1
    error('Input image should be normalized between 0 and 1.')
end

a   = 0.055;
thr = 0.0031308;  % threshold in LINEAR domain (not 0.04045)

I_enc = zeros(size(I_lin));

% Linear segment (low luminance)
mask = I_lin <= thr;
I_enc(mask)  = 12.92 * I_lin(mask);

% Power-law segment
I_enc(~mask) = (1 + a) * I_lin(~mask).^(1/2.4) - a;

end