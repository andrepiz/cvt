function [phi1, phi2, hphi, err] = sample_sphere_projected_uniform(alpha, nhphi, phimax)
% This function sample points between the boundaries of an
% illuminated sphere at phase angle alpha such that the projection of the
% arcs spanned by the angle intervals is constant on the projection plane.
% If phimax is not specified, pi/2 is used as maximum angle.

if ~exist('phimax','var')
    phimax = pi/2;
end

% Procedure: 
% solve the equation linking the consecutive angles with the projection 
% of the underlined spherical arc on the perpendicular plane 

% 2 * R * sin((phi2-phi1)/2) * cos(alpha - phi2 + (phi2-phi1)/2) = C

% Parametrizing phi2 = phi1 + hphi

% 2 * R * sin(hphi/2) * cos(alpha - phi1 - hphi/2) = C

% Solving iteratively:
% 1. phi1 = phi1min
% 2. find hphi
% 3. phi2 = phi1 + hphi
% 4. start again from 2.

% When phi1 = phi1min and phi2 = phi2max, C is equal to the projected illuminated horizon 
% C = 2 * R * sin(phi_span/2) * cos(alpha - phi2max + phi_span/2)

% So if we want to discretize in n_sectors, we will use C/n
% 2 * R * sin(hphi/2) * cos(alpha - phi1 - hphi/2) = 1/n * (2 * R * sin(phi_span/2) * cos(alpha - phi2max + phi_span/2))

% Dividing by 2 * R and re-arranging:
% n * (sin(hphi/2) * cos(alpha - phi1 - hphi/2)) - sin(phi_span/2) * cos(alpha - phi2max + phi_span/2) = 0

if alpha >= 0
    phi1min = max(alpha - pi/2, -phimax);
    phi2max = phimax;
else
    phi1min = -phimax;
    phi2max = min(alpha + pi/2, phimax);
end

phi_span = phi2max - phi1min;

% Objective functions
fun_obj_phi1 = @(hphi, alpha, phi1, n) n*sin(hphi/2).*cos(alpha - phi1 - hphi/2) - sin(phi_span/2) * cos(alpha - phi2max + phi_span/2);
fun_obj_phi2 = @(hphi, alpha, phi2, n) n*sin(hphi/2).*cos(alpha - phi2 + hphi/2) - sin(phi_span/2) * cos(alpha - phi2max + phi_span/2);

% init
phi1 = zeros(1, nhphi);
phi2 = zeros(1, nhphi);
hphi = zeros(1, nhphi);

phi1(1) = phi1min;
phi2(end) = phi2max;
hphi_temp = phi_span/nhphi;

if alpha >= 0
    % Forward
    for i = 1:nhphi-1
        fun_zero = @(x) fun_obj_phi1(x, alpha, phi1(i), nhphi);
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
    err = abs(sum(hphi)) - abs(phi_span);

elseif alpha < 0
    % Backward
    for i = nhphi:-1:2
        fun_zero = @(x) fun_obj_phi2(x, alpha, phi2(i), nhphi);
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
    err = abs(sum(hphi)) - abs(phi_span);

end

end