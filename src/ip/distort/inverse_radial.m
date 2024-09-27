function [UV, UV_c] = inverse_radial(RC, principalPoint, pxSize, radialDistortion)
cx = principalPoint(1); %[mm]
cy = principalPoint(2); %[mm]
f = 1; %cameraParams.K(1,1)*pxSize;

% RC are in image reference frame (origin in top-left corner, x row, y col)
RC_c = RC - [cx, cy];   % centered in optical center

% Transform in world units (mm) and normalize
XY = RC_c*pxSize/f;

% Define coefficients
k1 = radialDistortion(1); %in normalized image coordinates, world units
k2 = radialDistortion(2);

a1 = -k1;
a2 = 3*k1^2-k2;
a3 = 8*k1*k2-12*k1^3;
a4 = 5*(11*k1^4-11*k2*k1^2+k2^2);
a5 = -13*k1*(21*k1^4-28*k2*k1^2+6*k2^2);
a = [a1, a2, a3, a4, a5];
n = 5;

XYu_c = zeros(size(XY));
for i = 1:size(XY,1)
    x = XY(i,2);
    y = XY(i,1);
    r2 = x^2+y^2;
    % alpha = atan2(y,x);
    % p = [r2, r2^2, r2^3, r2^4, r2^5];
    % r_u = sqrt(r2)*(1+sum(a(1:n).*p(1:n)));
    % x_u = r_u*cos(alpha);
    % y_u = r_u*sin(alpha);

    x_u = x+x*(k1*r2+k2*r2^2);
    y_u = y+y*(k1*r2+k2*r2^2);
    XYu_c(i,:) = [y_u, x_u]; % [mm], centered in optical center
end

UV_c = XYu_c./pxSize*f; % in px, centered in optical center
UV = UV_c + [cx, cx];  % in px, image coordinates

% UV_2 = undistortPoints(XY, cameraParams)+0*[s/2+0.5 s/2+0.5];
end


 

