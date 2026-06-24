function [ang_tangency] = find_spheres_tangent_angle(d, R1, R2)
%FIND_SPHERES_TANGENT_ANGLE Given the tangent between two
%spheres of radius R1 and R2, find the angle from the line linking their
%centers where the tangent touches the sphere of radius R2. This angle is
%larger than 90 deg when R1 is bigger than R2, and lower than 90 deg when
%R1 is smaller than R2.

ang_inc = find_inclination_tangent(R1, R2, d);
ang_tangency = pi/2 + ang_inc;


end

