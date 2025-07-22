function img = map2grayscale(map, nbit, domain)
% Function to convert a map to a grayscale image. Latitude decreases with 
% increasing rows and longitude increases with
% increasing columns
% Input:
%   map [lat x lon] matrix where each element contains a value.
% Output:
%   rgb [u x v] grayscale image

% Scale the values into [0, 1]
ana = (map - domain(1)) ./ (domain(2) - domain(1));

img = analog2digital(ana, 1, 1, nbit);

end

