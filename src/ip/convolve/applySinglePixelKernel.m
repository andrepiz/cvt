function matout = applySinglePixelKernel(u, v, matin, kernel, flag_kernel_padding)
% This function applies a kernel to a specific pixel (u, v) in the matrix.
% The result is the convolution of the kernel at the pixel location.
% If flag_kernel_padding is true, padd the kernel with zeros up to the matin
% size before convolving.
%
% Inputs:
%   u, v   - the coordinates of the pixel where the kernel is applied
%   matin  - the input matrix
%   kernel - the kernel to be applied
%
% Output:
%   imgout - the filtered matrix

matout = matin;

if isempty(kernel)
    return
end

% Get the sizes
[lmat, hmat] = size(matin);
[lkern, hkern] = size(kernel);

usizekern = (lkern-1)/2;
vsizekern = (hkern-1)/2;

uminkern = u - usizekern;
umaxkern = u + usizekern;
vminkern = v - vsizekern;
vmaxkern = v + vsizekern;

if flag_kernel_padding
    kernel = padarray(kernel, [hmat, lmat], 'both');
    matout = applySinglePixelKernel(u, v, matin, kernel, false);
    return
end

rows_mat = max(1, uminkern):min(lmat, umaxkern); 
cols_mat = max(1, vminkern):min(hmat, vmaxkern);  

rows_shift = 1 - uminkern;
cols_shift = 1 - vminkern;

rows_kern = rows_mat + rows_shift;
cols_kern = cols_mat + cols_shift;

matout(rows_mat, cols_mat) = matin(rows_mat, cols_mat).*kernel(rows_kern, cols_kern);

end