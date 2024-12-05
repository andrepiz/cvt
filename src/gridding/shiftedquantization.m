function [bins, counts, edges] = shiftedquantization(coords, values, limits, granularity, params)
% Group and sum scattered data points into fixed-size quantiles with
% specified granularity.
% Each point value is kept in the nearest-neighbour sector or
% interpolated with a different method
% If shift is specified, perform the mean of the
% quantiles when the points are shifted of the specified size in all
% directions. If weight is specified, the mean is weighted with a kernel of
% size [NxN] where N is the maximum shift

arguments
    coords double
    values double
    limits double = [floor(min(coords,[],2)), 1 + ceil(max(coords,[],2))]
    granularity (1,1) double = 1
    params.method char = 'gaussian'
    params.flag_progress (1,1) logical = false
    params.shift (1,1) double = 0
    params.weight double = 1;
end

% Preliminary Checks
sz = (2*params.shift + 1);
nw = length(params.weight(:));
if nw == 1
    weight = params.weight*ones(sz);
    nw = length(weight(:));
else
    weight = params.weight;
end
if sz^2 ~= nw
    error('The weighting kernel must be of size [shift + 1 x shift + 1]')
end

weight_normalized = weight/sum(weight(:));
[centerrow, centercol] = ind2sub(sz, (nw + 1)/2);

bins_shift = cell(1, nw);
counts_shift = cell(1, nw);
for ix = 1:nw
    [ixrow, ixcol] = ind2sub(sz, ix);
    rowshift = ixrow - centerrow;
    colshift = ixcol - centercol;
    coords_shifted = coords + [rowshift; colshift];
    [bins_temp, counts_shift{ix}] = quantization(coords_shifted, values, limits, granularity, params);
    bins_shift{ix} = weight_normalized(ix)*bins_temp;
end

bins = sum(cat(3, bins_shift{:}), 3);
counts = sum(cat(3, counts_shift{:}), 3);
edges = cell(1, 2);
shift = min(0, limits(:, 1) - 1);
for kk = 1:2
    edges{kk} = [0:1:size(bins, kk)] + shift(kk)*granularity;
end


end

