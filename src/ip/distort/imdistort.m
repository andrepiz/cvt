function Idist = imdistort(I, principalPoint, pxSize, radialDistortionCoeff, UV)
%--------------------------------------------------------------------------
% Idist = distortImage(I, cameraParams, pxSize, radialDistortion) function
% that applies radial distortion to an input image I
%--------------------------------------------------------------------------
% Inputs:
%  I                      Undistorted image
%
%  cameraParams           A cameraParameters object obtained with the syntax 
%                          cameraParams = cameraParameters(paramStruct)
%
%  pxSize                 Pixel pitch [mm]
% 
%  radialDistortionCoeff       A 2-element vector [k1 k2]. 
%                          
%  UV                      Inverse mapping from input meshgrid, computed as
%                          [UV, ~] = inverse_radial(RC,cameraParams, pxSize, radialDistortion);
%--------------------------------------------------------------------------
% Output:
%  Idist                  Distorted image
%--------------------------------------------------------------------------

if nargin==4
% Create meshgrid from image 
nrows = size(I,1);
ncols = size(I,2);
[xi, yi] = meshgrid(1:ncols, 1:nrows);
RC = reshape([xi yi],[],2);

% Compute inverse mapping from input meshgrid
fprintf('\nCreating inverse mapping...')
tic
[UV, ~] = radial_distortion(RC, principalPoint, pxSize, radialDistortionCoeff);
tdist = toc;
fprintf('\nInverse mapping created in %fs',tdist)
end 

% Image distortion
fprintf('\nImage warping...')
tic
ifcn = @(c) UV;
tform = geometricTransform2d(ifcn);
Idist = imwarp(I,tform);
timwarp = toc;
fprintf('\nImage warped in %fs',timwarp)
end
