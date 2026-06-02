function img = encodeImage(img, encoding)
%ENCODEIMAGE Encode image given the encoding type

switch encoding
    case 'linear'
        return
        % do nothing

    case 'srgb'
        img = lin2srgb(img);

    otherwise
        error('Encoding not supported. Use linear or srgb.')
end

end