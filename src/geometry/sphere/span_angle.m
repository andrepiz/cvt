function [angle, c] = span_angle(dmin, dmax, delta)
% span_angle computes the span angle and the length of side c in a triangle
% given:
%   dmin   - length of side a (scalar or array)
%   dmax   - length of side b (scalar or array)
%   delta  - additional angle added to 90 degrees (in radians)
%
% Outputs:
%   angle  - computed angle between sides a and b (in radians)
%   c      - computed length of side c opposite this angle

% Total angle beta between sides a and c (90 degrees + delta)
beta = delta + pi/2;

a = dmin;
b = dmax;

% Compute side c using the derived formula:
% c = a*cos(beta) + sqrt(b^2 - a^2*sin^2(beta))
c = a.*cos(beta) + sqrt(b.^2 - (a.^2).*(sin(beta)).^2);

% Compute angle between sides a and b using law of cosines:
% angle = acos((a^2 + b^2 - c^2) / (2ab))
angle = acos((a.^2 + b.^2 - c.^2) ./ (2 .* a .* b));
end
