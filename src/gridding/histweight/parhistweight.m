function [bins, counts, edges] = parhistweight(coords, values, limits, granularity, method, threads)
% PARHISTWEIGHT parallel computing version of HISTWEIGHT with specified
% number of threads.
% HISTWEIGHT weights and bin scattered data points into uniform quantiles of 
% specified granularity within the specified limits. 
% Each data point is expressed in D-dimensional coordinates and has an 
% associated 1-dimensional value. The limits can be different for each dimension. 
% The granularity downsample the limits and increase the number of quantiles. 
% Each value is spread to the neighbouring 
% bins with a weight defined by three different methods:
%   invsquared: inverse squared distance with each vertex
%   diff: 1 minus distance normalized over maximum distance
%   area: fraction of square box area going to each sector
%
% INPUTS:
%   coords [D x N]
%   values [1 x N]
%   limits [D x 2]
%   granularity [1]
%   method [char]
%   threads [1]
%
% OUTPUTS:
%   bins [M1 x ... x Mi x ... MD]
%   counts [M1 x ... x Mi x ... MD]
%   edges [M1 x ... x Mi x ... MD]

% Preliminary checks
[D, N] = size(coords);
% intvs = round(linspace(0, N, threads+1));
% idxs_left = intvs(1:end-1) + 1;
% idxs_right = intvs(2:end);
spl = floor(N/threads);
rem = N - spl*threads;
if rem == 0
    rem = [];
end
coords_split = mat2cell(coords, D, [spl*ones(1, threads), rem]);
values_split = mat2cell(values, 1, [spl*ones(1, threads), rem]);

parfor ix = 1:length(coords_split)
    [bins_array(:,:,ix), counts_array(:,:,ix), edges(:,:,ix)] = histweight(coords_split{ix}, values_split{ix}, limits, granularity, method);
end

bins = sum(bins_array, 3);
counts = sum(counts_array, 3);
edges = edges(:,:,1);

end

