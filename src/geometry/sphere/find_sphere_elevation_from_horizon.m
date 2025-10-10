function el = find_sphere_elevation_from_horizon(d, R, phase_angle)  
% FIND_SPHERE_ELEVATION_FROM_HORIZON Computes the elevation at equator 
% from the horizon level of a point placed at distance d 
% and given phase angle from a sphere of radius R

Rd = R./d;
el = asin((sin(abs(phase_angle)) - Rd)./sqrt(1 + Rd.^2 -2.*Rd.*sin(abs(phase_angle)))); 

end