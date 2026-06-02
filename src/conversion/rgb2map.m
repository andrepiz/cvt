function img = rgb2map(img, nbit, encoding)
% Function to convert a rgb image to a normal map.
% Input:
%   img - HxWx3 RGB image with domain [0 2^nbit-1]
% Output:
%   img - HxWx3 matrix where each pixel contains the [nx, ny, nz] normal
%         vector with domain [-1 1]

if ~exist('encoding','var')
    encoding = 'linear';
end

img = digital2analog(img, 1, nbit, 1);
img = decodeImage(img, encoding);
img = img*2 - 1; % Scale the normal vectors from [0, 1] to [-1, 1]

end

