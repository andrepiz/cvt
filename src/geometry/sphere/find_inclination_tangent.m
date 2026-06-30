function ang_inc = find_inclination_tangent(r1, r2, d)
%FIND_INCLINATION_TANGENT Angle of the external tangent line between two spheres 
% with respect to the horizontal line connecting their two centers.
%
% Inputs:
%   r1 : radius of first sphere
%   r2 : radius of second sphere
%   d  : distance between centers
%
% Output:
%   ang_inc : angle [rad] of the external tangent line wrt the horizontal
%
% Notes:
%   For external tangent: sin(alpha) = (r1 - r2) / d
%   For internal tangent: sin(alpha) = (r1 + r2) / d

if d <= 0
    error('Distance d must be positive.');
end
if r1 < 0 || r2 < 0
    error('Radii must be non-negative.');
end
if abs(r1 - r2) > d
    error('Spheres are contained within each other; no external tangent exists.');
end

ang_inc = asin((r1 - r2) / d);

end