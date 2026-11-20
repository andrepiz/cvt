function noise = noise_fpn(ec_in, ec_fpn, pattern)
% Fixed-Pattern Noise enters the detector as an offset on the readout with 
% three possible patterns:
% - row: alternating offset values between rows
% - col: same as row, but on columns
% - checkerboard: like a checkerboard

sz = size(ec_in);
noise = ec_fpn(2)*ones(size(ec_in));
if size(ec_fpn, 1) > 1 && sz(2) > 1
    pattern = 'image';
    if size(ec_fpn, 1) == sz(1) && size(ec_fpn, 2) == sz(2)
    else
        error('When FPN is provided as a matrix, it should match the size of the image')
    end
elseif length(ec_fpn) < 2
    error('FPN should be provided as a vector of at least two elements or as a matrix of the same size of the image')
end

switch pattern
    case 'row'
        ixsA = 1:2:sz(2);
        noise(ixsA, :) = ec_fpn(1);

    case {'col','column'}
        ixsA = 1:2:sz(2);
        noise(:, ixsA) = ec_fpn(1);

    case 'checkerboard'
        ixsA = 1:2:sz(2);
        noise(ixsA, ixsA) = ec_fpn(1);
        noise(ixsA+1, ixsA+1) = ec_fpn(1);

    case 'image'
        noise = ec_fpn;

    otherwise
        error('Pattern not recognized. Use row, col or checkerboard');
end

end

