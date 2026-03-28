function [psf, psf_resized] = psf_zernike(wavelength, aberr_coeffs, pupil_diameter, focal_length, res_px, pixel_pitch, Nfft)
% COMPUTE_SINGLE_PSF  Fourier-optics PSF physically scaled to detector pixels.
%
% Inputs:
%   wavelength      - Single wavelength [m]
%   aberr_coeffs    - [1x9] Zernike coefficients [waves RMS], OSA/ANSI j=1..9
%   pupil_diameter  - Aperture diameter [m]
%   focal_length    - Focal length [m]
%   pixel_pitch     - Detector pixel pitch [m]
%   Nfft            - FFT grid size
%
% Output:
%   psf  - PSF kernel resampled to detector pixel scale, energy-normalized


    % The number of elements in the pupil function array is the same as the number of elements in the image plane array.
    % For an even number of array elements, the spacing between discrete coordinate values is equal to fS / numPixels
    % fS is called the sampling frequency and is equal to one divided by the spacing between image space coordinates. 
    % 
    fS = 1 ./ max(pixel_pitch); % Spatial sampling frequency [1/m]
    df = fS ./ max(res_px); % Spacing between discrete frequency coordinates [1/(m * px)]
    % fx = -fS / 2:df:fS / 2-df; % Spatial frequency [1/m]
    fNum = focal_length / pupil_diameter;
    fNA           = fNum / wavelength; % radius of the pupil [1/m]
    pupilSize     = fNA / df;        % [px]
    
    scale = Nfft / pupilSize;

    % --- Pupil grid (normalized coords) ---
    x_p = pupilSize/2*linspace(-1, 1, Nfft);  % [m]
    [xp, yp] = meshgrid(x_p);
    rho   = sqrt(xp.^2 + yp.^2)/(pupilSize/2);
    theta = atan2(yp, xp);
    aperture = double(rho <= 1);
    
    % --- Zernike wavefront ---
    coeffs = zeros(1, 9);
    if ~isempty(aberr_coeffs)
        n = min(numel(aberr_coeffs), 9);
        coeffs(1:n) = aberr_coeffs(1:n);
    end
    
    Z = zeros(Nfft, Nfft, 9);
    Z(:,:,1) = rho .* cos(theta);
    Z(:,:,2) = rho .* sin(theta);
    Z(:,:,3) = 2*rho.^2 - 1;
    Z(:,:,4) = rho.^2 .* sin(2*theta);
    Z(:,:,5) = rho.^2 .* cos(2*theta);
    Z(:,:,6) = (3*rho.^3 - 2*rho) .* sin(theta);
    Z(:,:,7) = (3*rho.^3 - 2*rho) .* cos(theta);
    Z(:,:,8) = rho.^3 .* sin(3*theta);
    Z(:,:,9) = rho.^3 .* cos(3*theta);
    
    phase_aberr = 2 * pi * sum(reshape(coeffs, 1, 1, 9) .* Z, 3);
    pupil       = aperture .* exp(1i * phase_aberr);
    
    psf_complex = fftshift(fft2(ifftshift(pupil))) * df^2;
    psf_full     = abs(psf_complex).^2;
    
    % Bilinear interpolation makes strict non-negative the PSF
    psf_resized = imresize(psf_full, ceil(res_px*scale/2)*2 + 1, 'bilinear');

    % Reduce further convolution by keeping ensquared energy only
    psf_ensq = crop_or_pad_psf(psf_resized, 0.99);
    
    % Normalize PSF
    psf = psf_ensq / sum(psf_ensq(:));
end

% -------------------------------------------------------------------------
function out = crop_or_pad_psf(psf, threshold_ensquared_energy)
% Keep only the central region with meaningful energy (99.9%), odd-sized.
    total    = sum(psf(:));
    target   = threshold_ensquared_energy * total;
    cx       = floor(size(psf,2)/2) + 1;
    cy       = floor(size(psf,1)/2) + 1;

    % Grow half-width until enclosing target energy
    for hw = 1:min(cx,cy)-1
        patch = psf(cy-hw:cy+hw, cx-hw:cx+hw);
        if sum(patch(:)) >= target
            break;
        end
    end
    out = psf(cy-hw:cy+hw, cx-hw:cx+hw);
end