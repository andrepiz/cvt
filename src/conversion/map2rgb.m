function img = map2rgb(img, nbit, encoding)
% Function to convert a 3D map to an RGB image. Latitude decreases with 
% increasing rows and longitude increases with increasing columns
% Input:
%   img - HxWx3 matrix where each pixel contains the [nx, ny, nz] normal
%         vector with domain [-1 1]
% Output:
%   img - HxWx3 RGB image with domain [0 2^nbit-1]

if ~exist('encoding','var')
    encoding = 'linear';
end

img = (img + 1) / 2; % Scale the normal vectors from [-1, 1] to [0, 1]
img = encodeImage(img, encoding);
img = analog2digital(img, 1, 1, nbit);

end

