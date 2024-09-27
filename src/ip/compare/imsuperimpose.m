function fh = imsuperimpose(background, varargin)

fh = figure();
grid on, hold on

imshow(background)
hold on
c = colorbar;
clim([0, max(background(:))])
xlabel('u [px]')
ylabel('v [px]')
c.Label.String = 'DN';

for ix = 1:length(varargin)
    img_binarized = varargin{ix};
    img_binarized(img_binarized > 1) = 255;
    img_binarized_rgb = zeros(size(img_binarized,1), size(img_binarized,2), 3);
    img_binarized_rgb(:,:,ix) = img_binarized;
    im_temp = imshow(img_binarized_rgb);
    im_temp.AlphaData = 0.3;
%    lg{ix} = ['Image', num2str(ix)];
end

%legend(lg)
title('Mask of image(s) superimposed on background image');

end