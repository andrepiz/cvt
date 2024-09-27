function img_psf = apply_psf(img, psf)
% img_psf = apply_psf(psf) function that apply a PSF to an image
% Inputs:
%   img     [UxV] Image to process
%   psf     [MxM] PSF to apply
%--------------------------------------------------------------------------
% Output:
%   img_psf [UxV] Processed image
%--------------------------------------------------------------------------

% resize, apply psf, resize
scale = size(psf, 1);
A = imresize(img, scale);
B = imfilter(A, psf);
img_psf = imresize(B, 1/scale);

end
