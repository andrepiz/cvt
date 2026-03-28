function tiffCorrect(filepath_tiff, filepath_tiff_new, data, bit_depth, scaling_factor)

% Create copy of file
if ~isfolder(fileparts(filepath_tiff_new))
    mkdir(fileparts(filepath_tiff_new))
end
copyfile(filepath_tiff, filepath_tiff_new);

% Define TIFF object
try
t = Tiff(filepath_tiff_new, 'w');
catch
t = Tiff(filepath_tiff_new, 'w8');
end

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
if exist("scaling_factor","var")
    if isscalar(scaling_factor)
        t.setTag('ImageDescription', sprintf('scaling=%d', scaling_factor));  % Standard metadata string
    elseif length(scaling_factor) == 2
        t.setTag('ImageDescription', sprintf('min=%f, max=%d', scaling_factor(1), scaling_factor(2)));  % Standard metadata string
    else
        error('Scaling factor should be a one-element or two-element (domain) vector')
    end 
end

% Write image and close file
t.write(data);
t.close();

end