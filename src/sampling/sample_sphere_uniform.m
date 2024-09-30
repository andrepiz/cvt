function [phi1, phi2, hphi] = sample_sphere_uniform(alpha, nhphi, phimax)
% This function sample points between the boundaries of an
% illuminated sphere at phase angle alpha with uniform sample step.

if ~exist('phimax','var')
    phimax = pi/2;
end

if alpha >= 0
    phi1min = max(alpha - pi/2, -phimax);
    phi2max = phimax;
else
    phi1min = -phimax;
    phi2max = min(alpha + pi/2, phimax);
end

hphi = (phi2max - phi1min)/(nhphi + 1);
phi = phi1min:hphi:phi2max;

phi1 = phi(1:end-1);
phi2 = phi(2:end);

end