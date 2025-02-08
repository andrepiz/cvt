function [phi1, phi2, hphi] = sample_sphere_uniform(phase_angle, nhphi, philims)
% This function sample points between the boundaries of an
% illuminated sphere at given phase angle with uniform sample step.
% If philims are not specified, [-pi/2 pi/2] are used as limits.

if ~exist('philims','var')
    philims = pi/2*[-1 1];
end
phimin = philims(1);
phimax = philims(2);

if phase_angle >= 0
    phi1min = max(phase_angle - pi/2, phimin);
    phi2max = phimax;
else
    phi1min = phimin;
    phi2max = min(phase_angle + pi/2, phimax);
end

hphi = (phi2max - phi1min)/(nhphi);
phi = phi1min:hphi:phi2max;

phi1 = phi(1:end-1);
phi2 = phi(2:end);

end