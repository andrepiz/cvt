function [phi1, phi2, hphi] = sample_sphere_uniform(phase_angle, nhphi, philims)
% This function samples points between the boundaries of an
% illuminated sphere at given phase angle with uniform sample step.
% If philims are not specified, [-pi/2 pi/2] are used as limits.

if ~exist('philims','var')
    philims = pi/2*[-1 1];
end
phi1min = philims(1);
phi2max = philims(2);

hphi = (phi2max - phi1min)/(nhphi);
phi = phi1min:hphi:phi2max;

phi1 = phi(1:end-1);
phi2 = phi(2:end);

end