function [P, u] = intersect_line_sphere(P1, P2, PC, R)
% Find the point of intersection of a line defined by two points P1 and P2
% with a sphere of radius R whose center is at PC. If two points of
% intersection are found, the one closer to P1 is returned.

P1 = reshape(P1, 3, []);
P2 = reshape(P2, 3, []);
PC = reshape(PC, 3, []);

a = vecnorm(P2 - P1).^2;
b = 2*((P2(1, :) - P1(1, :)).*(P1(1, :) - PC(1, :)) + ...
       (P2(2, :) - P1(2, :)).*(P1(2, :) - PC(2, :)) + ...
       (P2(3, :) - P1(3, :)).*(P1(3, :) - PC(3, :)));
c = vecnorm(PC).^2 + vecnorm(P1).^2 - ...
    2*(PC(1, :).*P1(1, :) + PC(2, :).*P1(2, :) + PC(3, :).*P1(3, :)) - ...
    R^2;

k = b.^2 - 4*a*c;

% the closest to P1
u = min((-b + sqrt(k))./(2*a), (-b - sqrt(k))./(2*a));
u(k < 0) = nan;
u(k == 0) = -b/(2*a);

P = P1 + u.*(P2 - P1);

end