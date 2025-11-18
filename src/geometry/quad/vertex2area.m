function quad_area = vertex2area(A, B, C, D)

% A, B, C, D: 3xN arrays, each column is a vertex in 3D
% quad_area: 1xN vector of quadrilateral areas

% Vectors for first triangle (A, B, C)
AB = B - A;
AC = C - A;

% Cross product and area for first triangle
cross1 = cross(AB, AC, 1);
area1 = 0.5 * sqrt(sum(cross1.^2, 1));

% Vectors for second triangle (A, C, D)
AD = D - A;

% Cross product and area for second triangle
cross2 = cross(AC, AD, 1);
area2 = 0.5 * sqrt(sum(cross2.^2, 1));

% Sum areas to get quadrilateral area
quad_area = area1 + area2;

end
