function [phi1, phi2, hphi, err] = sample_sphere_projected_uniform(phase_angle, nhphi, philims)
% This function sample points between the boundaries of an
% illuminated sphere at given phase angle such that the projection of the
% arcs spanned by the angle intervals is constant on the projection plane.
% If philims are not specified, [-pi/2 pi/2] are used as limits.

if ~exist('philims','var')
    philims = pi/2*[-1 1];
end
phimin = philims(1);
phimax = philims(2);

% Procedure: 
% solve the equation linking the consecutive angles with the projection 
% of the underlined spherical arc on the perpendicular plane 

% 2 * R * sin((phi2-phi1)/2) * cos(phase_angle - phi2 + (phi2-phi1)/2) = C

% Parametrizing phi2 = phi1 + hphi

% 2 * R * sin(hphi/2) * cos(phase_angle - phi1 - hphi/2) = C

% Solving iteratively:
% 1. phi1 = phi1min
% 2. find hphi
% 3. phi2 = phi1 + hphi
% 4. start again from 2.

% When phi1 = phi1min and phi2 = phi2max, C is equal to the projected illuminated horizon 
% C = 2 * R * sin(phi_span/2) * cos(phase_angle - phi2max + phi_span/2)

% So if we want to discretize in n_sectors, we will use C/n
% 2 * R * sin(hphi/2) * cos(phase_angle - phi1 - hphi/2) = 1/n * (2 * R * sin(phi_span/2) * cos(phase_angle - phi2max + phi_span/2))

% Dividing by 2 * R and re-arranging:
% n * (sin(hphi/2) * cos(phase_angle - phi1 - hphi/2)) - sin(phi_span/2) * cos(phase_angle - phi2max + phi_span/2) = 0

% Find the number of points to sample with projected-uniform sampling
% We find left and right limits due to phase angle visibility.
% What is missing is sampled with uniform sampling
phi_span = phimax - phimin;
phi1min = max(phase_angle - pi/2, phimin);
phi1rem = min(0, phimin - (phase_angle - pi/2));
n1rem = ceil(abs(phi1rem)/phi_span*nhphi);
n1rem(isnan(n1rem)) = 0;
phi2max = min(phase_angle + pi/2, phimax);
phi2rem = max(0, phimax - (phase_angle + pi/2));
n2rem = ceil(phi2rem/phi_span*nhphi);
n2rem(isnan(n2rem)) = 0;
nhphi_red = nhphi - n1rem - n2rem;
phi_span_red = phi2max - phi1min;

% Objective functions
fun_obj_phi1 = @(hphi, phase_angle, phi1, n) n*sin(hphi/2).*cos(phase_angle - phi1 - hphi/2) - sin(phi_span_red/2) * cos(phase_angle - phi2max + phi_span_red/2);
fun_obj_phi2 = @(hphi, phase_angle, phi2, n) n*sin(hphi/2).*cos(phase_angle - phi2 + hphi/2) - sin(phi_span_red/2) * cos(phase_angle - phi2max + phi_span_red/2);
%fun_obj_phi1 = @(hphi, phase_angle, phi1, n) n*sin(hphi/2).*cos(phi1 + hphi/2 - phase_angle) - sin(phi_span_red/2) * cos(phi2max - phi_span_red/2 - phase_angle);
%fun_obj_phi2 = @(hphi, phase_angle, phi2, n) n*sin(hphi/2).*cos(phi2 - hphi/2 - phase_angle) - sin(phi_span_red/2) * cos(phi2max - phi_span_red/2 - phase_angle);

% init
phi1 = zeros(1, nhphi_red);
phi2 = zeros(1, nhphi_red);
hphi = zeros(1, nhphi_red);

phi1(1) = phi1min;
phi2(end) = phi2max;
hphi_temp = phi_span_red/nhphi_red;

if phase_angle >= 0
    % Forward
    for i = 1:nhphi_red-1
        fun_zero = @(x) fun_obj_phi1(x, phase_angle, phi1(i), nhphi_red);
        [hphi(i), ~] = fzero(fun_zero, hphi_temp);
        if hphi(i) < 0
            % correction in case step is negative
            hphi(i) = -hphi(i);
        end
        %rad2deg(hphi(i))
        phi2(i) = phi1(i) + hphi(i);
        phi1(i+1) = phi2(i);
        hphi_temp = hphi(i);
    end
    hphi(end) = phi2(end) - phi1(end);
    err = abs(sum(hphi)) - abs(phi_span_red);

elseif phase_angle < 0
    % Backward
    for i = nhphi_red:-1:2
        fun_zero = @(x) fun_obj_phi2(x, phase_angle, phi2(i), nhphi_red);
        [hphi(i), ~] = fzero(fun_zero, hphi_temp);
        if hphi(i) < 0
            % correction in case step is negative
            hphi(i) = -hphi(i);
        end
        %rad2deg(hphi(i))
        phi1(i) = phi2(i) - hphi(i);
        phi2(i-1) = phi1(i);
        hphi_temp = hphi(i);
    end
    hphi(1) = phi2(1) - phi1(1);
    err = abs(sum(hphi)) - abs(phi_span_red);

end

% Fill with uniform sampling the points outside the limits
[phi1left, phi2left, hphileft] = sample_sphere_uniform(0, n1rem, [phimin phi1min]);
[phi1right, phi2right, hphiright] = sample_sphere_uniform(0, n2rem, [phi2max phimax]);

% Create final vector
phi1 = [phi1left, phi1, phi1right];
phi2 = [phi2left, phi2, phi2right];
hphi = [hphileft, hphi, hphiright];

end