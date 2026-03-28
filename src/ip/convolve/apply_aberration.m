function [ec_out, psf] = apply_aberration(ec_in, pupil_diameter, focal_length, pixel_pitch, wavelength, aberr_coeffs, varargin)
% APPLY_ABERRATION  End-to-end optical aberrations with physically-scaled diffraction PSF.
%
% Inputs (required):
%   ec_in          - [N x M]     electron count matrix, panchromatic [electrons]
%                  - [N x M x W] per-channel, W must match numel(wavelength)
%   pupil_diameter - scalar    aperture diameter [m]
%   focal_length   - scalar    focal length [m]
%   pixel_pitch    - scalar    detector pixel pitch [m]
%   wavelength     - [W x 1]   wavelength(s) [m], scalar = mono, vector = polychromatic
%   aberr_coeffs   - [1x9]     Zernike coefficients [waves RMS], OSA/ANSI j=1..9
%                              [Z1_tiltX, Z2_tiltY, Z3_defocus, Z4_astig45,
%                               Z5_astig0, Z6_comaY, Z7_comaX, Z8_trefoilY, Z9_trefoilX]
%                              Pass [] for diffraction-limited (all zeros)
%
% Inputs (optional Name-Value):
%   'ChannelWeights'  - [W x 1] spectral weights per wavelength [ones(W,1)]
%                               Automatically normalized to sum = 1.
%   'PadMode'         - conv2 padding: 'same' | 'full' | 'valid'    ['same']
%   'Oversample'      - sub-pixel oversampling factor (integer)      [1]
%   'Nfft'            - Size of PSF FTT       [1]
%
% Output:
%   ec_out  - [N x M]     blurred (panchromatic input)
%           - [N x M x W] per-channel blurred (multichannel input)
%
% Notes:
%   - For polychromatic light, PSFs are weighted by 'Weights'.
%     If not specified, all wavelengths contribute equally.
%
% Example (polychromatic, uniform weights):
%   ec_blur = apply_aberration(ec_sharp, 0.05, 1.0, 5e-6,          ...
%               [500e-9, 550e-9, 600e-9, 650e-9, 700e-9], aberr);
%
% Example (polychromatic, custom weights):
%   w = [0.1 0.3 0.5 0.3 0.1];
%   ec_blur = apply_aberration(ec_sharp, 0.05, 1.0, 5e-6,          ...
%               [500e-9, 550e-9, 600e-9, 650e-9, 700e-9], aberr,   ...
%               'Weights', w);
%
% Author: Andrea Pizzetti — March 2026

    wavelength = wavelength(:);   % ensure column vector
    W = numel(wavelength);
    
    % ADD: detect multichannel input
    is_multichannel = (ndims(ec_in) > 1);
    if is_multichannel
        assert(size(ec_in, 3) == W, ...
            'The number of color channels of the image (%d) must be equal to the number of wavelengths (%d).', size(ec_in,3), W);
    end
    
    % ADD: validate and expand aberr_coeffs to [W x 9]
    if isempty(aberr_coeffs)
        aberr_coeffs = zeros(W, 9);
    elseif size(aberr_coeffs, 1) == 1
        aberr_coeffs = repmat(aberr_coeffs(:)', W, 1);
    else
        assert(size(aberr_coeffs,1) == W && size(aberr_coeffs,2) == 9, ...
            'Aberration coefficients must be [1x9] or [Wx9] where W is the number of wavelengths (%d).', W);
    end

    % --- Input parsing ---
    p = inputParser;
    addParameter(p, 'ChannelWeights', ones(numel(wavelength), 1), ...
            @(x) isvector(x) && numel(x) == numel(wavelength) && all(x >= 0));
    addParameter(p, 'PadMode',       'same',              @ischar);
    addParameter(p, 'Oversample',     1,                  @(x) isscalar(x) && x >= 1);
    addParameter(p, 'Nfft',           512,                  @isscalar);

    parse(p, varargin{:});

    channel_weights = p.Results.ChannelWeights;
    pad_mode        = p.Results.PadMode;
    os_factor       = round(p.Results.Oversample);
    Nfft            = p.Results.Nfft;

    % --- Oversampling ---
    if os_factor > 1
        ec_in_os = repelem(ec_in, os_factor, os_factor);
    else
        ec_in_os = ec_in;
    end
    if is_multichannel && os_factor > 1
        ec_in_os = repelem(ec_in, os_factor, os_factor, 1);
    end

    res_px = size(ec_in, [1 2]);

    % --- PSF computation ---
    if is_multichannel

        ec_out_os = zeros(size(ec_in_os));
        psf       = cell(W, 1);   % store per-channel PSF for output
        parfor i = 1:W
            psf{i}           = psf_zernike(wavelength(i), aberr_coeffs(i,:), ...
                                           pupil_diameter, focal_length, res_px, pixel_pitch, Nfft);
            ec_out_os(:,:,i) = conv2(ec_in_os(:,:,i), psf{i}, pad_mode);
        end

    else

        if isscalar(wavelength)
        
            psf = psf_zernike(wavelength, aberr_coeffs(1,:), ...         
                              pupil_diameter, focal_length, res_px, pixel_pitch, Nfft);
        
        else

            weights = channel_weights / sum(channel_weights);
            parfor i = 1:W
                psf{i}           = psf_zernike(wavelength(i), aberr_coeffs(i,:), ...
                                           pupil_diameter, focal_length, res_px, pixel_pitch, Nfft);
            end
            sz1       = max(cellfun(@(x) size(x, 1), psf));
            sz2       = max(cellfun(@(x) size(x, 2), psf));
            psf_accum = zeros(sz1, sz2);
            for i = 1:W
                psf_accum = match_and_add(psf_accum, psf{i}, weights(i), [sz1, sz2]);
            end
            psf = psf_accum / sum(psf_accum(:));
        end

        % --- Convolution ---
        ec_out_os = conv2(ec_in_os, psf, pad_mode);
    end


    if os_factor > 1
        ec_out = imresize(ec_out_os, 1/os_factor, 'bilinear');
    else
        ec_out = ec_out_os;
    end
    
    if is_multichannel && os_factor > 1
        ec_out = zeros(size(ec_in));
        parfor i = 1:W
            ec_out(:,:,i) = imresize(ec_out_os(:,:,i), 1/os_factor, 'bilinear');
        end
    end

end

% -------------------------------------------------------------------------
function out = match_and_add(accum, psf_new, w, target_sz)
% Pad both arrays to target_sz before accumulating, to handle size mismatches.
    accum_pad = pad_to_size(accum,   target_sz);
    new_pad   = pad_to_size(psf_new, target_sz);
    out       = accum_pad + w * new_pad;
end

function out = pad_to_size(psf, target_sz)
% Zero-pad PSF symmetrically to target_sz (must be >= size(psf)).
    dp  = target_sz - size(psf);
    dp  = max(dp, 0);
    out = padarray(psf, floor(dp/2), 0, 'pre');
    out = padarray(out, ceil(dp/2),  0, 'post');
end