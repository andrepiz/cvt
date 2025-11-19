function [phi1, phi2, hphi] = sample_sphere_concentrated(phase_angle, nhphi, philims, concentrationFactor)
% This function samples points between the boundaries of a
% sphere with sample points more concentrated around the nadir point.
%
% Inputs:
%   nhphi: number of sample intervals
%   philims: 2-element vector with [min_phi max_phi] limits (default [-pi/2 pi/2])
%
% Outputs:
%   phi1, phi2: interval edges of sampled phi
%   hphi: average sample step

if ~exist('philims','var') || isempty(philims)
    philims = pi/2*[-1 1];
end

if ~exist('concentrationFactor','var')
    concentrationFactor = 5;
end

phi_min = philims(1);
phi_max = philims(2);
phi_ref = max(phi_min, min(phi_max, phase_angle));

% Normalize limits and reference to [0,1]
phiNorm_min = 0;
phiNorm_max = 1;
phiNorm_ref = (phi_ref - phi_min) / (phi_max - phi_min);

% Create a uniform grid in [0,1]
u = linspace(phiNorm_min, phiNorm_max, nhphi+1);

A = 1 / (1 + exp(concentrationFactor * phiNorm_ref));
B = 1 / (1 + exp(-concentrationFactor * (1 - phiNorm_ref)));

numerator = 1 - u * (B - A) - A;
denominator = u * (B - A) + A;

% To avoid numerical issues clamp inside log argument
clip_limit = 1e-12;
ratio = max(clip_limit, numerator ./ denominator);

y = phiNorm_ref - (1 / concentrationFactor) * log(ratio);

phi = phi_min + y*(phi_max - phi_min);

% Map back to original phi range

phi1 = phi(1:end-1);
phi2 = phi(2:end);

hphi = phi2 - phi1;

end
