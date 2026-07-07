function out = poissrnd_fast(in)
%POISSRND_APPROX Fast Poisson distribution generator approximating to a
%gaussian for large Poisson rates

out = in;
mask = in > 20;
out(~mask) = poissrnd(in(~mask));
% For larger poisson rate >10/20, poission distribution converges to gaussian
out(mask)  = in(mask) + sqrt(in(mask)).*randn(nnz(mask),1);
out = max(round(out), 0);

end

