function psf = create_gaussian_psf(sigma, support)
% psf = create_gaussian_psf(sigma, support) create a Gaussian PSF with a
% specified standard deviation and support (size).
% Inputs:
%   sigma     [1] 
%   support   [2] 
%--------------------------------------------------------------------------
% Output:
%   psf       [MxN] 
%--------------------------------------------------------------------------

[x, y] = meshgrid(-floor(support(1)/2):floor(support(1)/2), -floor(support(2)/2):floor(support(2)/2));
psf = (1/(2*pi*sigma^2)) * exp(-(x.^2 + y.^2) / (2*sigma^2));


end
