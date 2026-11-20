%% Spatial LUMIO-Cam PSFs: nine illuminated squares
% Run this script from any directory. Results remain in the workspace:
% img, img_psf, psfs, metrics, elapsed_seconds, relative_flux_error.
% Angles preserve the Zemax signs, with positive image directions right/down
% as required by apply_spatialPsf. No detector noise or blooming is added.
% Width and skewness describe the combined optical response; they do not
% isolate Zernike coma from other aberrations or from finite pixel sampling.

cvt_install()
psf_folder = fullfile(cvt_home, 'test', 'data', 'psf');

%% Load the 450 nm, 3-by-3 grid (columns: x; rows: y)
node_tokens = {'m0.0262', 'p0.0000', 'p0.0262'};
kernel_filenames = cell(3);
for v = 1:3
    for u = 1:3
        kernel_filenames{v,u} = fullfile(psf_folder, sprintf( ...
            ['LUMIO-Cam_PSFkernel_wl_450.0000nm_bxoh_%s_byoh_%s_' ...
             'fp_VIS_samp_2.166666um_relative.txt'], node_tokens{u}, node_tokens{v}));
    end
end
% Use the known design sampling rather than the rounded 2.167 um header.
[psfs, spacing, offset_horizontal, offset_vertical] = ...
    load_spatialPsf(kernel_filenames, 2.167e-6, [], [], 'zemax');

%% Build nine isolated, equally bright 3-by-3 pixel squares
image_size = [1024, 1024];
muPixel = 13e-6;            % Detector pixel pitch [m]
focal_length = 0.200;      % Test geometry [m]; keeps all nine targets in frame
square_width = 9;          % Odd integer; set to 1 to test point sources
square_intensity = 1;      % Integrated signal per illuminated pixel
assert(mod(square_width,2)==1 && square_width>=1);
half_square = (square_width-1)/2;
image_centre = (image_size+1)/2;
target_cols = round(image_centre(2) + focal_length*tan(offset_horizontal)/muPixel);
target_rows = round(image_centre(1) + focal_length*tan(offset_vertical)/muPixel);

native_sizes = cellfun(@(p) max(size(p)), psfs);
half_window = half_square + ceil(max(native_sizes(:))*spacing/muPixel) + 3;
assert(all(target_cols-half_window>=1) && all(target_cols+half_window<=image_size(2)) ...
    && all(target_rows-half_window>=1) && all(target_rows+half_window<=image_size(1)), ...
    'Targets and their measurement windows must fit inside the image.');
assert(all(abs(diff(target_cols))>2*half_window) ...
    && all(abs(diff(target_rows))>2*half_window), 'Measurement windows must not overlap.');
img = zeros(image_size);
for v = 1:3
    for u = 1:3
        img(target_rows(v)+(-half_square:half_square), ...
            target_cols(u)+(-half_square:half_square)) = square_intensity;
    end
end

%% Apply the spatial PSFs (pool startup is excluded from the reported time)
if license('test','Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
    parpool;
end
timer = tic;
img_psf = apply_spatialPsf(img, psfs, offset_horizontal, offset_vertical, ...
    focal_length, spacing, muPixel);
elapsed_seconds = toc(timer);
% First-call worker initialization/data transfer can still affect this timing.
relative_flux_error = abs(sum(img_psf(:))-sum(img(:)))/sum(img(:));
assert(isequal(size(img_psf),size(img)) && all(isfinite(img_psf(:))));
assert(min(img_psf(:))>=-1e-12, 'PSF produced negative signal.');
assert(relative_flux_error<1e-10, 'Signal was lost despite isolated, interior targets.');

%% Measure each square in an isolated window
% Coordinates are relative to the geometric centre of each input square.
% Added RMS removes the known square's second moment; skewness uses central
% third moments, so a simple image translation is not mistaken for skewness.
expected_flux = square_width^2*square_intensity;
square_variance = (square_width^2-1)/12;
[xx, yy] = meshgrid(-half_window:half_window);
measurements = zeros(9,11);
spot_patches = cell(3);
for v = 1:3
    for u = 1:3
        k = (v-1)*3+u;
        patch = img_psf(target_rows(v)+(-half_window:half_window), ...
                        target_cols(u)+(-half_window:half_window));
        spot_patches{v,u} = patch;
        flux = sum(patch(:));
        probability = patch/flux;
        dx = sum(probability(:).*xx(:));
        dy = sum(probability(:).*yy(:));
        cx = xx-dx; cy = yy-dy;
        variance_x = sum(probability(:).*cx(:).^2);
        variance_y = sum(probability(:).*cy(:).^2);
        added_rms = sqrt(max(0,variance_x+variance_y-2*square_variance));
        skew_x = sum(probability(:).*cx(:).^3)/max(variance_x,eps)^(3/2);
        skew_y = sum(probability(:).*cy(:).^3)/max(variance_y,eps)^(3/2);
        % This source-centred asymmetry includes centroid shift as well.
        inverted = rot90(patch,2);
        asymmetry = sum(abs(patch(:)-inverted(:)))/(2*flux);
        actual_x = atan((target_cols(u)-image_centre(2))*muPixel/focal_length);
        actual_y = atan((target_rows(v)-image_centre(1))*muPixel/focal_length);
        measurements(k,:) = [v,u,rad2deg(actual_x),rad2deg(actual_y), ...
            flux/expected_flux,dx,dy,added_rms,skew_x,skew_y,asymmetry];
    end
end
metrics = array2table(measurements, 'VariableNames', ...
    {'GridRow','GridCol','FieldX_deg','FieldY_deg','FluxFraction', ...
     'ShiftX_px','ShiftY_px','AddedRMS_px','SkewX','SkewY','Asymmetry'});
assert(all(abs(metrics.FluxFraction-1)<1e-10), 'A measurement window missed signal.');
fprintf('Spatial PSF application: %.3f s (pool startup excluded).\n',elapsed_seconds);
fprintf('Relative total flux error: %.3g\n',relative_flux_error);
disp(metrics);
fprintf(['Shift: centroid displacement. AddedRMS: blur after removing square size.\n' ...
    'SkewX/Y: directional tail asymmetry about the centroid.\n' ...
    'Asymmetry: source-centred inversion difference, including displacement.\n' ...
    'These diagnose coma-like effects; they do not identify individual aberration modes.\n']);

%% Full scene: equal intensity scales and signed residual
figure('Name','Spatial PSF - nine squares','Color','w');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; imagesc(img,[0 square_intensity]); axis image; colorbar;
title('Input: nine equal squares'); xlabel('Column'); ylabel('Row');
nexttile; imagesc(img_psf,[0 square_intensity]); axis image; colorbar;
title(sprintf('Spatial PSF: %.3f s',elapsed_seconds)); xlabel('Column'); ylabel('Row');
nexttile;
residual = img_psf-img;
limit = max(max(abs(residual(:))),eps);
imagesc(residual,[-limit limit]); axis image; colorbar;
title('Output - input'); xlabel('Column'); ylabel('Row');
colormap(parula);

%% Enlarged detector responses: common log scale exposes faint tails
figure('Name','Spatial PSF - square responses','Color','w');
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for v = 1:3
    for u = 1:3
        k = (v-1)*3+u;
        nexttile;
        imagesc(-half_window:half_window,-half_window:half_window, ...
            log10(max(spot_patches{v,u}/square_intensity,1e-6)),[-6 0]);
        axis image; hold on;
        rectangle('Position',[-square_width/2,-square_width/2,square_width,square_width], ...
            'EdgeColor','w','LineStyle',':');
        plot(0,0,'w+');
        plot(metrics.ShiftX_px(k),metrics.ShiftY_px(k),'rx','LineWidth',1.5);
        title(sprintf('(%+.1f,%+.1f) deg | RMS %.3f px\nskew (%+.2f,%+.2f)', ...
            rad2deg(offset_horizontal(u)),rad2deg(offset_vertical(v)), ...
            metrics.AddedRMS_px(k),metrics.SkewX(k),metrics.SkewY(k)));
        xlabel('Offset x [pixels]'); ylabel('Offset y [pixels]'); colorbar;
    end
end
sgtitle('log10(signal / input intensity); white box: source, red x: centroid');

%% Native PSFs: inspect optical structure before detector integration
figure('Name','Spatial PSF - native Zemax kernels','Color','w');
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for v = 1:3
    for u = 1:3
        nexttile;
        kernel = psfs{v,u};
        x_um = ((1:size(kernel,2))-(floor(size(kernel,2)/2)+1))*spacing*1e6;
        y_um = ((1:size(kernel,1))-(floor(size(kernel,1)/2)+1))*spacing*1e6;
        imagesc(x_um,y_um,log10(max(kernel/max(kernel(:)),1e-6)),[-6 0]);
        axis image; hold on; plot(0,0,'w+');
        title(sprintf('(%+.1f,%+.1f) deg', ...
            rad2deg(offset_horizontal(u)),rad2deg(offset_vertical(v))));
        xlabel('Offset x [um]'); ylabel('Offset y [um]'); colorbar;
    end
end
sgtitle('Native PSFs: log10(intensity / local peak), before pixel integration');
