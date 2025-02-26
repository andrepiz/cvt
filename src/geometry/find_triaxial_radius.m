function R = find_triaxial_radius(lon, lat, R3)
% FIND_TRIAXIAL_RADIUS calculates the radius at a given latitude and longitude
% for a sphere, biaxial ellipsoid, or triaxial ellipsoid.
%
% Inputs:
%   lon    - Longitude [rad]
%   lat    - Latitude [rad]
%   R3     - Tri-Axial Ellipsoid Radii. 
%            Semi-major axis (equatorial, x-direction)
%            Intermediate axis (equatorial, y-direction) [optional]
%            Semi-minor axis (polar, z-direction) [optional]
%
% Output:
%   R      - Radius at the given latitude and longitude

if size(lat, 1) ~= size(lon, 1) || size(lat, 2) ~= size(lon, 2)
    error('Please provide same dimensions for latitude and longitude coordinates')
end

switch length(R3)
    case 1 % Sphere (only mean radius given)
        R = R3*ones(size(lat));
    
    case 2 % Biaxial Ellipsoid (assume b = a, c is given)
        a = R3(1);
        c = R3(2);
        R = (a * c) ./ sqrt((c^2 .* cos(lat).^2) + (a^2 .* sin(lat).^2));
    
    case 3 % Triaxial Ellipsoid (full model with a, b, and c)
        a = R3(1);
        b = R3(2);
        c = R3(3);
        R = (a * b * c) ./ sqrt((b^2 .* c^2 .* cos(lat).^2 .* cos(lon).^2) + ...
                                (a^2 .* c^2 .* cos(lat).^2 .* sin(lon).^2) + ...
                                (a^2 .* b^2 .* sin(lat).^2));
    
    otherwise
        error('Please provide radius as a vector of 1 to 3 components');
end

end