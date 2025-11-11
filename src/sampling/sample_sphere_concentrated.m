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

phi_min = philims(1);
phi_max = philims(2);
phi_ref = max(phi_min, min(phi_max, phase_angle));

% Normalize limits and reference to [0,1]
phiNorm_min = 0;
phiNorm_max = 1;
phiNorm_ref = (phi_ref - phi_min) / (phi_max - phi_min);

% Create a uniform grid in [0,1]
u = linspace(phiNorm_min, phiNorm_max, nhphi+1);

% method = 'exponential';
% switch lower(method)
% 
%     case 'exponential'
%         lambda = 10; % concentration strength parameter >0
%         warp_func = @(x) phiNorm_ref + sign(x - phiNorm_ref) .* (1/lambda) .* (1 - exp(-lambda*abs(x - phiNorm_ref)));
% 
%         v = warp_func(u);
%         vNorm = (v - v(1))./(v(end) - v(1));
% 
%     case 'logistic'

        A = 1 / (1 + exp(concentrationFactor * phiNorm_ref));
        B = 1 / (1 + exp(-concentrationFactor * (1 - phiNorm_ref)));
        
        numerator = 1 - u * (B - A) - A;
        denominator = u * (B - A) + A;
        
        % To avoid numerical issues clamp inside log argument
        clip_limit = 1e-12;
        ratio = max(clip_limit, numerator ./ denominator);
        
        y = phiNorm_ref - (1 / concentrationFactor) * log(ratio);

        phi = phi_min + y*(phi_max - phi_min);

%     case 'quadratic'
%         beta = 0.3; % concentration strength, positive number
%         warp_func = @(x) max(0, min(1, x + beta*((x - phiNorm_ref).^2).*(x < phiNorm_ref).* (-1) + beta*((x - phiNorm_ref).^2).*(x > phiNorm_ref)));
% 
%         vNorm = warp_func(u);
% 
%     case 'gaussian'
%         sigma = 0.1; % width of Gaussian bump
%         gamma = 0.15; % magnitude
%         gauss = @(x) gamma * exp(-((x - phiNorm_ref).^2) / (2*sigma^2));
%         warp_func = @(x) max(0, min(1, x + gauss(x)));
% 
%     otherwise
%         error('Unknown warping method');
% end


% % Ensure monotonic increasing (clip any deviation)
% v = sort(v);

% Map back to original phi range

phi1 = phi(1:end-1);
phi2 = phi(2:end);

hphi = phi2 - phi1;

end
