function kern = gaussianKernel(kern_size, kern_sigma, flag_plot)

[xk, yk] = meshgrid(-kern_size:kern_size, -kern_size:kern_size);
arg = -(xk.*xk + yk.*yk) ./ (2 * kern_sigma.*kern_sigma);
kern = exp(arg); 
kern = kern / max(kern(:));

if flag_plot
    figure()
    imagesc(kern)
    colorbar
end

end