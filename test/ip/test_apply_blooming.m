%% BLOOMING Test  Visualize sphere rendering with varying blooming params
cvt_install()
clear; close all; clc;

cm = colormap('turbo');

img_size = 512*[1 1];

%% PSF
fNum = 1.4;
focal_length = 0.025;
pupil_diameter = focal_length/fNum; 
pixel_pitch = 2e-6;
wavelength = 600e-9;
aberr = zeros(1,9); aberr(3) = 2; 

% PSF
[psf_ensquared, psf] = psf_zernike(wavelength, aberr, pupil_diameter, focal_length, img_size, pixel_pitch, 1);

figure()
surf(psf,'EdgeColor','none'); set(gca(), 'ZScale','log')

%% 1. SATURATION
[ecr_ideal] = render_test_annulus(img_size);

fwc = 100e3;
tExp_for_sat = fwc/max(ecr_ideal(:));
tExp = tExp_for_sat*70;

ec_in = ecr_ideal*tExp;
[ec_out, trapped, psf] = apply_blooming(ec_in, fwc, 'gaussian');
%[ec_out, trapped, psf] = apply_blooming(ec_in, fwc, 'gaussian_2d');
%[ec_out, trapped, psf] = apply_blooming(ec_in, fwc, psf);
ec_out_noblooming = min(ec_in, fwc);

figure()
subplot(1,2,1)
imagesc(ec_out_noblooming)
axis image
colormap('gray')
colorbar
subplot(1,2,2)
imagesc(ec_out)
axis image
colormap('gray')
colorbar

figure()
imagesc(imabsdiff(ec_out_noblooming, ec_out))
axis image
colormap(cm)
colorbar

%%
function ec = render_test_annulus(sz)
% RENDER_TEST_ANNULUS  Generate a monochrome annulus (disk with hole).
%
% Inputs:
%   sz  [1x2]  image size [rows, cols]
%
% Output:
%   ec  [rows x cols]  monochrome image scaled to electron counts (0-1000)

    [X, Y] = meshgrid(1:sz(2), 1:sz(1));

    cx = sz(2) / 2;
    cy = sz(1) / 2;

    r_outer = min(sz) / 3;
    r_inner = min(sz) / 6;   % hole radius

    dist = sqrt((X - cx).^2 + (Y - cy).^2);

    mask = (dist <= r_outer) & (dist >= r_inner);

    ec = 1000 * double(mask);
end