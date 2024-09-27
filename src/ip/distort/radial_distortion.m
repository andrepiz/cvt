function [UV, UV_c] = radial_distortion(RC, principalPoint, pxSize, radialDistortionCoeff)
cx = principalPoint(1); %[px]
cy = principalPoint(2); %[px]
f = 1; 

% RC are in image reference frame (origin in top-left corner, x row, y col)
RC_c = RC - [cx, cy];   % [px] centered in optical center

% Transform in world units (mm) and normalize
XY = RC_c*pxSize/f; %[mm]

% Define coefficients
k1 = radialDistortionCoeff(1); %in normalized image coordinates, world units
k2 = radialDistortionCoeff(2);


XYu_c = zeros(size(XY));
for i = 1:size(XY,1)
    x = XY(i,1);
    y = XY(i,2);
    r2 = x^2+y^2;

    x_u = x+x*(k1*r2+k2*r2^2);
    y_u = y+y*(k1*r2+k2*r2^2);
    XYu_c(i,:) = [x_u, y_u]; % [mm], centered in optical center
end

UV_c = XYu_c./pxSize*f; % [px] centered in optical center
UV = UV_c + [cx, cx];  % [px] image coordinates

end

