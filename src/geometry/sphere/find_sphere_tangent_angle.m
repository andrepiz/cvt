function [ang_tangency, ang_bearing] = find_sphere_tangent_angle(d, R)
%FIND_SPHERE_TANGENT_ANGLE find the angle of tangency for an observer at
%distance d looking at a sphere with radius R. If the observer is inside
%the sphere, return nan.

ang_bearing = asin(R./d);
ang_bearing(d<R) = nan;
ang_tangency = pi/2 - ang_bearing;

end

