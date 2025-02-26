function tiffResize(filepath_tiff, filepath_tiff_resized, resize_factor)

data = imread(filepath_tiff);

data = imresize(data, resize_factor);

tiffCorrect(filepath_tiff, filepath_tiff_resized, data, whos(data).bytes)

end