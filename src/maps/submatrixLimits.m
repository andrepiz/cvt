%======================================================================
% submatrixLimits
%   Return row‑ and column‑index limits for an N×N tiling of a matrix.
%
% SYNTAX
%   [rowStart,rowEnd,colStart,colEnd] = submatrixLimits(M, K, N)
%   [rowStart,rowEnd,colStart,colEnd] = submatrixLimits(A,     N)
%
% INPUT
%   M, K : (optional) dimensions of the big matrix.
%   A    : (alternative) the matrix itself.
%   N    : positive integer – number of blocks along each dimension.
%
% OUTPUT (all N×N double matrices)
%   rowStart(i,j) = first row   of tile (i,j)
%   rowEnd  (i,j) = last  row   of tile (i,j)
%   colStart(i,j) = first column of tile (i,j)
%   colEnd  (i,j) = last  column of tile (i,j)
%
%   Tile (i,j) therefore corresponds to rows rowStart(i,j):rowEnd(i,j)
%   and columns colStart(i,j):colEnd(i,j), using 1‑based MATLAB indexing.
%======================================================================
function [rowStart, rowEnd, colStart, colEnd] = submatrixLimits(varargin)

    % -------- 1. Parse inputs -------------------------------------------
    if nargin == 2                  % usage: (A, N)
        [M, K] = size(varargin{1});
        N      = varargin{2};
    elseif nargin == 3              % usage: (M, K, N)
        M = varargin{1};
        K = varargin{2};
        N = varargin{3};
    else
        error(['Usage: submatrixLimits(A, N)  or ', ...
               'submatrixLimits(M, K, N)']);
    end
    if ~(isscalar(N) && N == round(N) && N > 0)
        error('N must be a positive integer.');
    end

    % -------- 2. Determine strip sizes (rows & columns) ------------------
    rBase  = floor(M / N);
    rSizes = rBase + (1:N <= mod(M, N));    % first rem rows get +1

    cBase  = floor(K / N);
    cSizes = cBase + (1:N <= mod(K, N));    % first rem cols get +1

    rEdge = [0, cumsum(rSizes)];            % row boundaries, length N+1
    cEdge = [0, cumsum(cSizes)];            % col boundaries, length N+1

    % -------- 3. Build N×N matrices of limits ---------------------------
    rowStart = repmat(rEdge(1:N).' + 1, 1, N);   % broadcast down columns
    rowEnd   = repmat(rEdge(2:N+1).',   1, N);

    colStart = repmat(cEdge(1:N)  + 1,  N, 1);   % broadcast across rows
    colEnd   = repmat(cEdge(2:N+1),     N, 1);
end
