function img_out = imreconstruct(img_in, granularity, kern, filter, filter_size)
% Reconstruct an image by downsampling it with a specified filter. If kern is
% specified, apply a convolution kernel before resizing.
%
% "nearest","box"            Nearest-neighbor interpolation; the output pixel is assigned the value of the pixel that the point falls within. No other pixels are considered.
% "bilinear","triangle"      Bilinear interpolation; the output pixel value is a weighted average of pixels in the nearest 2-by-2 neighborhood
% "bicubic","cubic"          Bicubic interpolation; the output pixel value is a weighted average of pixels in the nearest 4-by-4 neighborhood. Note: Bicubic interpolation can produce pixel values outside the original range.
% "lanczos2"	             Lanczos-2 kernel
% "lanczos3"	             Lanczos-3 kernel
% "lanczos"                  Lanczos kernel. Specify width parameter
% "gaussian"                 Gaussian kernel. Specify width and sigma parameters
% "osculatory"               Osculatory rational interpolation

if isempty(kern)
    img_conv = img_in;
else
    kern_energy_preserving = kern/(sum(kern(:)));
    img_conv = conv2(img_in, kern_energy_preserving, 'same');
end

if ~exist('filter','var')
    filter = 'lanczos3';
end

switch filter

    case 'gaussian'
            img_out = granularity^2*imresize(img_conv, 1/granularity, {@(x)gaussianResampling(x, filter_size/6), filter_size});

    case 'lanczos'
            img_out = granularity^2*imresize(img_conv, 1/granularity, {@(x) lanczosResampling(x, filter_size), filter_size});

    case 'osculatory'
            img_out = granularity^2*imresize(img_conv, 1/granularity, {@oscResampling, filter_size});

    otherwise
        
        try
            img_out = granularity^2*imresize(img_conv, 1/granularity, filter);
        catch
            error('Reconstruction filter not recognized')
        end

end

function f = gaussianResampling(x,a)
    absx = abs(x);
    absx2 = absx.^2;
    sigma = 1/(2*a);  % https://dsp.stackexchange.com/questions/75899/appropriate-gaussian-filter-parameters-when-resizing-image
    f = exp(-absx2./sigma^2);
    f(abs(x) > 3*sigma) = 0;
    f(x == 0) = 1;
end

function f = lanczosResampling(x,a)
    f = a*sin(pi*x) .* sin(pi*x/a) ./ ...
        (pi^2 * x.^2);
    f(abs(x) > a) = 0;
    f(x == 0) = 1;
end

function f = oscResampling(x)
    absx = abs(x);
    absx2 = absx.^2;
    
    f = (absx <= 1) .* ...
        ((-0.168*absx2 - 0.9129*absx + 1.0808) ./ ...
        (absx2 - 0.8319*absx + 1.0808)) ...
        + ...
        ((1 < absx) & (absx <= 2)) .* ...
        ((0.1953*absx2 - 0.5858*absx + 0.3905) ./ ...
        (absx2 - 2.4402*absx + 1.7676));
end

end