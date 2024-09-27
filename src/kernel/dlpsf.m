function [Irel] = dlpsf(r, lambda, fNum)
% The intensity of the Airy pattern follows the Fraunhofer diffraction pattern 
% of a circular aperture, given by the squared modulus of the Fourier transform 
% of the circular aperture

% I(theta) = I0*(2*J1(x)/x)^2
% Where J1(x) is the Bessel function of the first kind of order one
% given 
%   x = k R sin(theta)
%   k = 2*pi/lambda the wavenumber
%   R the radius of aperture
%   theta the angle of observation
%   
% Note: The Airy pattern is observable when f/lambda >> R^2 (i.e. in the far field)

% Parametrization with r = f*tan(theta) = f*sin(theta) for small theta
% Irel = I/I0 = (2*J1(xmod)/x)^2
% given
%   x = 2*pi/lambda*D/2*r/f = pi/lambda*r/fNum

x = pi/lambda*r/fNum;
Irel = (2*besselj(1, x)./x).^2;
end