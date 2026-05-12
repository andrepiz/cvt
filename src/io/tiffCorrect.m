function tiffCorrect(filepath_tiff, filepath_tiff_new, data, bit_depth, scaling_factor)

% Create copy of file
fld = fileparts(filepath_tiff_new);
if ~isfolder(fld) && ~isempty(fld)
    mkdir(fld)
end
if ~isequal(filepath_tiff, filepath_tiff_new)
    copyfile(filepath_tiff, filepath_tiff_new);
end

% Define TIFF object
try
    t = Tiff(filepath_tiff_new, 'w');
catch
    t = Tiff(filepath_tiff_new, 'w8');
end

% Set TIFF tags for floating-point storage
tagstruct.ImageLength = size(data, 1);
tagstruct.ImageWidth = size(data, 2);
nc = size(data, 3);
if nc == 1
    tagstruct.SamplesPerPixel = 1;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
elseif nc == 3
    tagstruct.SamplesPerPixel = 3;
    tagstruct.Photometric = Tiff.Photometric.RGB;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
else
    tagstruct.SamplesPerPixel = nc;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Separate;
    tagstruct.ExtraSamples      = repmat(Tiff.ExtraSamples.Unspecified, 1, nc-1);
end

tagstruct.BitsPerSample = bit_depth; % bit depth storage

if isa(data,'double') || isa(data, 'single')
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
elseif isa(data,'int8') || isa(data, 'int16') || isa(data, 'int32')
    tagstruct.SampleFormat = Tiff.SampleFormat.Int;
elseif isa(data,'uint8') || isa(data, 'uint16') || isa(data, 'uint32')
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
end

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