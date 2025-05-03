function tiffResize(filepath_tiff, filepath_tiff_resized, resize_factor)

data = imread(filepath_tiff);

data = imresize(data, resize_factor);

switch whos('data').class
    case 'uint8'
        nbytes = 8;
    case 'int16'
        nbytes = 16;
    case 'single'
        nbytes = 32;
    otherwise
        error('class not recognized')
end

tiffCorrect(filepath_tiff, filepath_tiff_resized, data, nbytes)

end