function [locus, cu, cv] = conic2locus(conic, np, cu, cv)
% Inputs:
%   conic - [A B C D E F] or 3x3 conic matrix
%   np - number of points in the locus

% Extract coefficients
if size(conic, 1) == 3 && size(conic, 2) == 3
    A = conic(1,1);
    B = 2*conic(1,2);
    C = conic(2,2);
    D = 2*conic(1,3);
    E = 2*conic(2,3);
    F = conic(3,3);
elseif length(conic) == 6
    A = conic(1); B = conic(2); C = conic(3);
    D = conic(4); E = conic(5); F = conic(6);
else
    error('Input a 3x3 matrix or 1x6 vector for conic');
end

if ~exist("cu", 'var') || ~exist("cv", 'var')
    % Calculate ellipse center (U0,V0)
    M = [2*A, B; B, 2*C];
    b = [-D; -E];
    center = M \ b;
    cu = center(1);
    cv = center(2);
end

% Translate origin to (U0,V0)
Fp = F + D*cu + E*cv + A*cu^2 + B*cu*cv + C*cv^2;

% Form quadratic form matrix Q
Q = [A, B/2; B/2, C];
[V_eig, D_eig] = eig(Q);
lambda = diag(D_eig);

% Handle based on type of conic (discriminant)
discriminant = B^2 - 4*A*C;

t = linspace(0, 2*pi, np)';

if discriminant < 0
    % Ellipse or circle
    % if Fp >= 0
    %     error('Fp must be negative for an ellipse');
    % end
    a = sqrt(-Fp / lambda(1));
    b = sqrt(-Fp / lambda(2));
    uv = V_eig * [a*cos(t)'; b*sin(t)'];
elseif discriminant > 0
    % Hyperbola
    s = linspace(-2, 2, np)';
    a = sqrt(abs(Fp / lambda(1)));
    b = sqrt(abs(Fp / lambda(2)));

    % Two branches
    uv1 = V_eig * [a*cosh(s)'; b*sinh(s)'];
    uv2 = V_eig * [-a*cosh(s)'; -b*sinh(s)'];
    uv = [uv1, uv2];  % Combine branches
else
    % Parabola: one eigenvalue ≈ 0
    % Sample values along u-axis and solve for v
    u = linspace(-5, 5, np);
    v = zeros(size(u));
    for i = 1:np
        % Solve quadratic in v: A u^2 + B u v + C v^2 + D u + E v + F = 0
        % At shifted origin (U0,V0), remove linear terms
        uu = u(i);
        coeffs = [C, B*uu, A*uu^2 + Fp];
        r = roots(coeffs);
        if ~isempty(r) && isreal(r)
            v(i) = r(1); % take one branch
        else
            v(i) = NaN;
        end
    end
    uv = [u; v];
end

% Translate back to original coordinates
U = uv(1,:) + cu;
V = uv(2,:) + cv;
locus = [U; V];

end
