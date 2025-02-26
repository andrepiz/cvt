function [LCS, L, C, S, Limg, Cimg, Simg] = ssim_comparison(img, tmpl, noise_level, flag_apply_ncc, flag_plot, dnr, cmap)

if ~exist('cmap','var')
    cmap = colormap('gray');
end

if ~exist('dnr','var')
    dnr = diff(getrangefromclass(tmpl));
end

% Remove unwanted pixels for comparison
imgmask_undernoise = img < noise_level;
tmplmask_undernoise = tmpl < noise_level;

if flag_apply_ncc

    % Compute Normalized Cross-Correlation on the images where we set equal the
    % background
    img_ncc = img;
    img_ncc(imgmask_undernoise) = 0;
    tmpl_ncc = tmpl;
    tmpl_ncc(tmplmask_undernoise) = 0;
    nccMatrix = normxcorr2(img_ncc, tmpl_ncc);
    
    % Find the peak of NCC
    [maxNccValue, maxIndex] = max(abs(nccMatrix(:)));
    [ypeak, xpeak] = ind2sub(size(nccMatrix), maxIndex);
    
    % Determine the offset
    offsetY = ypeak - size(img, 1);
    offsetX = xpeak - size(img, 2);
    
    % Align the image to the template
    img = imtranslate(img, [offsetX, offsetY],'FillValues',nan);
    imgmask_undernoise = imtranslate(imgmask_undernoise, [offsetX, offsetY]);
end


[L, Limg] = ssim(img, tmpl, 'exponents',[1 0 0],'DynamicRange', dnr,'Radius',0.5);
[C, Cimg] = ssim(img, tmpl, 'exponents',[0 1 0],'DynamicRange', dnr,'Radius',0.5);
[S, Simg] = ssim(img, tmpl, 'exponents',[0 0 1],'DynamicRange', dnr,'Radius',0.5);

% Remove unwanted pixels from SSIM comparison
mask_noise = imgmask_undernoise | tmplmask_undernoise;
Limg(mask_noise) = nan;
Cimg(mask_noise) = nan;
Simg(mask_noise) = nan;
L = mean(Limg(:), 'omitnan');
C = mean(Cimg(:), 'omitnan');
S = mean(Simg(:), 'omitnan');
LCS = L*C*S;

img_plot = double(img);
img_plot(imgmask_undernoise) = nan;
tmpl_plot = double(tmpl);
tmpl_plot(tmplmask_undernoise) = nan;

if flag_plot

    mask_img = img_plot/dnr;
    mask_img(mask_img>0) = uint8(255);
    mask_img(isnan(mask_img)) = 0;
    mask_img = repmat(mask_img, 1, 1, 3);
    mask_img(:,:, [2, 3]) = 0;

    mask_tmpl = tmpl_plot/dnr;
    mask_tmpl(mask_tmpl>0) = uint8(255);
    mask_tmpl(isnan(mask_tmpl)) = 0;
    mask_tmpl = repmat(mask_tmpl, 1, 1, 3);
    mask_tmpl(:,:, [1, 2]) = 0;
    
    mask_rgb = mask_img + mask_tmpl;

    figure('name','ncc_blue_template_mask_over_red_image_mask')
    h_mask = imshow(mask_rgb);
    h_mask.AlphaData = 0*(all(mask_rgb==0, 3)) + 1*(any(mask_rgb>0, 3));

    figure('Name','noise_mask_over_image')
    hold on
    mask_noise_plot = repmat(double(mask_noise), 1, 1, 3);
    mask_noise_plot(:, :, [2, 3]) = 0;
    h_nan = imshow(mask_noise_plot);
    h_nan.AlphaData = 0.5*any(mask_noise_plot, 3);
    hold on
    h_img = imshow(img_plot/max(img_plot,[],'all'));
    h_img.AlphaData = 0.5;

    clmin = min([min(Limg(:)),min(Cimg(:)),min(Simg(:))]);

    f1 = figure('name',['luminance_',num2str(L)]); 
    hL = imshow(Limg);
    colormap(f1, cmap)
    hL.AlphaData = 0*mask_noise + 1*~mask_noise;
    cb = colorbar();
    clim([clmin, 1])
    xlabel('u [px]')
    ylabel('v [px]')
    cb.Label.String = 'Similarity Index [-]';
    cb.TickLabelInterpreter = "latex";
    cb.Label.Interpreter = "latex";

    f2 = figure('name',['contrast_',num2str(C)]);
    hC = imshow(Cimg);
    colormap(f2, cmap)
    hC.AlphaData = 0*mask_noise + 1*~mask_noise;
    cb = colorbar();
    clim([clmin, 1])
    xlabel('u [px]')
    ylabel('v [px]')
    cb.Label.String = 'Similarity Index [-]';
    cb.TickLabelInterpreter = "latex";
    cb.Label.Interpreter = "latex";

    f3 = figure('name',['structure_',num2str(S)]);
    hS = imshow(Simg);
    colormap(f3, cmap)
    hS.AlphaData = 0*mask_noise + 1*~mask_noise;
    cb = colorbar();
    clim([clmin, 1])
    xlabel('u [px]')
    ylabel('v [px]')
    cb.Label.String = 'Similarity Index [-]';
    cb.TickLabelInterpreter = "latex";
    cb.Label.Interpreter = "latex";
end

end