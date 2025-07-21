function tiffCorrect(filepath_tiff, filepath_tiff_new, data, bit_depth)

% Create copy of file
if ~isfolder(fileparts(filepath_tiff_new))
    mkdir(fileparts(filepath_tiff_new))
end
copyfile(filepath_tiff, filepath_tiff_new);

% Define TIFF object
t = Tiff(filepath_tiff_new, 'w');

% Set TIFF tags for floating-point storage
tagstruct.ImageLength = size(data, 1);
tagstruct.ImageWidth = size(data, 2);
if size(data, 3) == 1
    tagstruct.SamplesPerPixel = 1;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
elseif size(data, 3) == 3
    tagstruct.SamplesPerPixel = 3;
    tagstruct.Photometric = Tiff.Photometric.RGB;
end
tagstruct.BitsPerSample = bit_depth; % bit depth storage
if isa(data,'double') || isa(data, 'single')
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
elseif isa(data,'int8') || isa(data, 'int16') || isa(data, 'int32')
    tagstruct.SampleFormat = Tiff.SampleFormat.Int;
elseif isa(data,'uint8') || isa(data, 'uint16') || isa(data, 'uint32')
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
end
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
t.setTag(tagstruct);

% Write image and close file
t.write(data);
t.close();

end