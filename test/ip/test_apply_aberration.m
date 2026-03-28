%% OPTICAL_BLURRING_TEST  Visualize sphere rendering with varying blur params
cvt_install()
clear; close all; clc;

cm = colormap('gray');

%% 1. ABERRATIONS 
img_size = 256*[1 1];
[ec_sharp] = render_test_annulus(img_size);

fNum = 1.4;
focal_length = 0.025;
pupil_diameter = focal_length/fNum; 
pixel_pitch = 2e-6;
wavelength = 600e-9;

% Diffraction-Limited
ec_diffraction_limited = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength, []);

% Tilt
aberr = zeros(1,9); 
aberr(1) = 15; aberr(2) = 19; % [Z1_tiltX, Z2_tiltY, Z3_defocus, Z4_astig45,
%                               Z5_astig0, Z6_comaY, Z7_comaX, Z8_trefoilY, Z9_trefoilX]
ec_tilt = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength, aberr);
fprintf(['\nTotal signal change: ', num2str((sum(ec_sharp(:)) - sum(ec_tilt(:)))./sum(ec_sharp(:))*1e2),'%%'])
fprintf(['\nRRMS signal change: ', num2str(rms(imabsdiff(ec_sharp, ec_tilt)./max(ec_sharp(:)), 'all')*1e2),'%%'])

% Blurring
aberr = zeros(1,9); aberr(3) = 3; 
ec_blurred = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength, aberr);
fprintf(['\nTotal signal change: ', num2str((sum(ec_sharp(:)) - sum(ec_blurred(:)))./sum(ec_sharp(:))*1e2),'%%'])
fprintf(['\nRRMS signal change: ', num2str(rms(imabsdiff(ec_sharp, ec_blurred)./max(ec_sharp(:)), 'all')*1e2),'%%'])

% Astigmatism
aberr = zeros(1,9); 
aberr(4) = 5; aberr(5) = 5; % [Z1_tiltX, Z2_tiltY, Z3_defocus, Z4_astig45,
%                               Z5_astig0, Z6_comaY, Z7_comaX, Z8_trefoilY, Z9_trefoilX]
ec_astigmatism = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength, aberr);
fprintf(['\nTotal signal change: ', num2str((sum(ec_sharp(:)) - sum(ec_astigmatism(:)))./sum(ec_sharp(:))*1e2),'%%'])
fprintf(['\nRRMS signal change: ', num2str(rms(imabsdiff(ec_sharp, ec_astigmatism)./max(ec_sharp(:)), 'all')*1e2),'%%'])

% Coma
aberr = zeros(1,9); 
aberr(6) = 8; aberr(7) = 8; % [Z1_tiltX, Z2_tiltY, Z3_defocus, Z4_astig45,
%                               Z5_astig0, Z6_comaY, Z7_comaX, Z8_trefoilY, Z9_trefoilX]
ec_coma = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength, aberr);
fprintf(['\nTotal signal change: ', num2str((sum(ec_sharp(:)) - sum(ec_coma(:)))./sum(ec_sharp(:))*1e2),'%%'])
fprintf(['\nRRMS signal change: ', num2str(rms(imabsdiff(ec_sharp, ec_coma)./max(ec_sharp(:)), 'all')*1e2),'%%'])

% Trifoil
aberr = zeros(1,9); 
aberr(8) = 5; aberr(9) = 5; % [Z1_tiltX, Z2_tiltY, Z3_defocus, Z4_astig45,
%                               Z5_astig0, Z6_comaY, Z7_comaX, Z8_trefoilY, Z9_trefoilX]
ec_trifoil = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength, aberr);
fprintf(['\nTotal signal change: ', num2str((sum(ec_sharp(:)) - sum(ec_trifoil(:)))./sum(ec_sharp(:))*1e2),'%%'])
fprintf(['\nRRMS signal change: ', num2str(rms(imabsdiff(ec_sharp, ec_trifoil)./max(ec_sharp(:)), 'all')*1e2),'%%'])

% All
aberr = zeros(1,9); 
aberr(1) = 15; aberr(2) = 19; 
aberr(3) = 3; aberr(4) = 5; 
aberr(5) = 5; aberr(6) = 8; 
aberr(7) = 8; aberr(8) = 5; 
aberr(9) = 5; 
ec_aberrated = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength, aberr);
fprintf(['\nTotal signal change: ', num2str((sum(ec_sharp(:)) - sum(ec_blurred(:)))./sum(ec_sharp(:))*1e2),'%%'])
fprintf(['\nRRMS signal change: ', num2str(rms(imabsdiff(ec_sharp, ec_blurred)./max(ec_sharp(:)), 'all')*1e2),'%%'])

figure('Name', 'Optical Blurring Test Suite', 'Position', [100 100 1700 900]);
subplot(2,4, 1); 
imagesc(ec_sharp); 
axis image
colormap(cm);
title('Ideal')

% Plot
subplot(2,4, 2);
imagesc(ec_diffraction_limited); 
axis image
colormap(cm);
cb = colorbar;
clm = cb.Limits;
clim(clm)
title('Diffraction-Limited')

% Plot
subplot(2,4, 3);
imagesc(ec_blurred); 
axis image
colormap(cm);
cb = colorbar;
clim(clm)
title('Blur')

% Plot
subplot(2,4, 4);
imagesc(ec_tilt); 
axis image
colormap(cm);
cb = colorbar;
clim(clm)
title('Tilt')

% Plot
subplot(2,4, 5);
imagesc(ec_astigmatism); 
axis image
colormap(cm);
cb = colorbar;
clim(clm)
title('Astigmatism')

% Plot
subplot(2,4, 6);
imagesc(ec_coma); 
axis image
colormap(cm);
cb = colorbar;
clim(clm)
title('Coma')

% Plot
subplot(2,4, 7);
imagesc(ec_trifoil); 
axis image
colormap(cm);
cb = colorbar;
clim(clm)
title('Trifoil')

% Plot
subplot(2,4, 8);
imagesc(ec_aberrated); 
axis image
colormap(cm);
% cb = colorbar;
clim(clm)
title('All')

%% 2. Nfft accuracy

aberr = zeros(1,9); aberr(3) = 3; 
ec_blurred_v1 = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength, aberr, 'Nfft', 512);
ec_blurred_v2 = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength, aberr, 'Nfft', 2048);
fprintf(['\nRRMS signal difference: ', num2str(rms(imabsdiff(ec_blurred_v1, ec_blurred_v2)./max(ec_blurred_v2(:)), 'all')*1e2),'%%'])

% Plot
figure()
subplot(1,2, 1);
imagesc(ec_blurred_v1); 
axis image
colormap(cm);
clim(clm)
subplot(1,2, 2);
imagesc(ec_blurred_v2); 
axis image
colormap(cm);
% cb = colorbar;
clim(clm)

%% 3. Wavelength

aberr = zeros(1,9); aberr(3) = 4; 
wavelength_1 = 100e-9;
ec_blurred_v1 = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength_1, aberr);
wavelength_2 = 900e-9;
ec_blurred_v2 = apply_aberration(ec_sharp, pupil_diameter, focal_length, pixel_pitch, wavelength_2, aberr);
fprintf(['\nRRMS signal difference: ', num2str(rms(imabsdiff(ec_blurred_v1, ec_blurred_v2)./max(ec_blurred_v2(:)), 'all')*1e2),'%%'])

% Plot
figure()
subplot(1,2, 1);
imagesc(ec_blurred_v1); 
axis image
colormap(cm);
clim(clm)
title(['Wavelength: ', num2str(wavelength_1*1e9),'nm'])
subplot(1,2, 2);
imagesc(ec_blurred_v2); 
axis image
colormap(cm);
title(['Wavelength: ', num2str(wavelength_2*1e9),'nm'])
% cb = colorbar;
clim(clm)

%% 4. fNum

aberr = zeros(1,9); aberr(3) = 4; 
fNum_1 = 1.4;
pupil_diameter_1 = focal_length/fNum_1; 
ec_blurred_v1 = apply_aberration(ec_sharp, pupil_diameter_1, focal_length, pixel_pitch, wavelength, aberr);
fNum_2 = 4;
pupil_diameter_2 = focal_length/fNum_2; 
ec_blurred_v2 = apply_aberration(ec_sharp, pupil_diameter_2, focal_length, pixel_pitch, wavelength, aberr);
fprintf(['\nRRMS signal difference: ', num2str(rms(imabsdiff(ec_blurred_v1, ec_blurred_v2)./max(ec_blurred_v2(:)), 'all')*1e2),'%%'])


% Plot
figure()
subplot(1,2, 1);
imagesc(ec_blurred_v1); 
axis image
colormap(cm);
clim(clm)
title(['f#: ', num2str(fNum_1),''])
subplot(1,2, 2);
imagesc(ec_blurred_v2); 
axis image
colormap(cm);
title(['f#: ', num2str(fNum_2),''])
% cb = colorbar;
clim(clm)

%%
function [ec_sharp, rgb_weights] = render_test_sphere(sz)
% RENDER_TEST_SPHERE  Generate a test image with three overlapping colored disks.
%
% Inputs:
%   sz   [1x2]  image size [rows, cols]
%
% Outputs:
%   ec_sharp    [rows x cols x 3]  RGB image scaled to electron counts (0-1000)
%   rgb_weights [1x3]              luminosity weights [R G B]

    [X, Y] = meshgrid(1:sz(2), 1:sz(1));

    % --- Disk centers and radius ---
    r  = min(sz) / 4;
    offset = r * 0.6;   % overlap amount

    cx1 = sz(2)/2;            cy1 = sz(1)/2 - offset;   % top disk    → Red
    cx2 = sz(2)/2 - offset;   cy2 = sz(1)/2 + offset;   % bottom-left → Green
    cx3 = sz(2)/2 + offset;   cy3 = sz(1)/2 + offset;   % bottom-right→ Blue

    % --- Binary disk masks ---
    mask1 = sqrt((X - cx1).^2 + (Y - cy1).^2) <= r;
    mask2 = sqrt((X - cx2).^2 + (Y - cy2).^2) <= r;
    mask3 = sqrt((X - cx3).^2 + (Y - cy3).^2) <= r;

    % --- Layered compositing (later disks paint over earlier ones) ---
    ec_sharp = zeros([sz(1) sz(2) 3]);

    % Disk 1: Red
    ec_sharp(:,:,1) = ec_sharp(:,:,1) + mask1 * 1.0;
    ec_sharp(:,:,2) = ec_sharp(:,:,2) + mask1 * 0.0;
    ec_sharp(:,:,3) = ec_sharp(:,:,3) + mask1 * 0.0;

    % Disk 2: Green (overwrites overlap with disk 1)
    ec_sharp(:,:,1) = ec_sharp(:,:,1) .* ~mask2 + mask2 * 0.0;
    ec_sharp(:,:,2) = ec_sharp(:,:,2) .* ~mask2 + mask2 * 1.0;
    ec_sharp(:,:,3) = ec_sharp(:,:,3) .* ~mask2 + mask2 * 0.0;

    % Disk 3: Blue (overwrites overlaps with disks 1 and 2)
    ec_sharp(:,:,1) = ec_sharp(:,:,1) .* ~mask3 + mask3 * 0.0;
    ec_sharp(:,:,2) = ec_sharp(:,:,2) .* ~mask3 + mask3 * 0.0;
    ec_sharp(:,:,3) = ec_sharp(:,:,3) .* ~mask3 + mask3 * 1.0;

    % --- Scale to electron counts ---
    ec_sharp = 1000 * ec_sharp;

    rgb_weights = [0.3 0.59 0.11];
end

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
