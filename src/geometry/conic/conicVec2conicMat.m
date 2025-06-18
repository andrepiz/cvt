function conicMat = conicVec2conicMat(conicVec)
% conicVec2conicMat builds a 3x3 conic matrix from [A, B, C, D, E, F]
%
% INPUT:
%   conic_vec - 1x6 or 6x1 vector: [A, B, C, D, E, F]
%
% OUTPUT:
%   conic_mat - 3x3 symmetric conic matrix

if numel(conicVec) == 5
    conicVec(6) = 1;
elseif numel(conicVec) ~= 6
    error('Input must be a 6-element vector: [A, B, C, D, E, F]');
end

A = conicVec(1);
B = conicVec(2);
C = conicVec(3);
D = conicVec(4);
E = conicVec(5);
F = conicVec(6);

conicMat = [ A,   B/2, D/2;
      B/2, C,  E/2;
      D/2, E/2, F  ];
end