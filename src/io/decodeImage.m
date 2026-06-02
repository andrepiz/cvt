function img = decodeImage(img, encoding)
%DECODEIMAGE Decode image given the encoding type

switch encoding
    case 'linear'
        return
        % do nothing

    case 'srgb'
        img = srgb2lin(img);

    otherwise
        error('Encoding not supported. Use linear or srgb.')
end

end