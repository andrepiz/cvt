function window = extractCenteredWindow(matrix, winrows, wincols)
% This function extracts a winrows x wincols window from the center of a matrix.
% If the window exceeds the matrix boundaries, it is padded with zeros.
%
% Inputs:
%   matrix  - the input matrix (MxN)
%   winrows - the number of rows of the window
%   wincols - the number of columns of the window
%
% Output:
%   window  - the extracted O x P window centered on the matrix

% Get the size of the input matrix
[M, N] = size(matrix);

% Compute the center of the matrix
centerRow = ceil(M / 2);
centerCol = ceil(N / 2);

% Determine the start and end indices for the rows and columns in the matrix
rowStart = centerRow - floor(winrows / 2);
rowEnd = centerRow + floor(winrows / 2);

colStart = centerCol - floor(wincols / 2);
colEnd = centerCol + floor(wincols / 2);

% Initialize the output window with zeros (O x P)
window = zeros(winrows, wincols);

% Calculate the valid row and column ranges that fall inside the matrix
rowStartValid = max(1, rowStart);
rowEndValid = min(M, rowEnd);

colStartValid = max(1, colStart);
colEndValid = min(N, colEnd);

% Determine the positions in the window where the matrix values will be placed
rowOffsetStart = rowStartValid - rowStart + 1;
colOffsetStart = colStartValid - colStart + 1;

% Copy the relevant portion of the matrix into the zero-padded window
window(rowOffsetStart:rowOffsetStart + (rowEndValid - rowStartValid), ...
       colOffsetStart:colOffsetStart + (colEndValid - colStartValid)) = ...
       matrix(rowStartValid:rowEndValid, colStartValid:colEndValid);

end