function kern = exponentialKernel(kern_size, kern_lambda, flag_plot)

[xk, yk] = meshgrid(-kern_size:kern_size, -kern_size:kern_size);
kern = exp(-kern_lambda*sqrt(xk.^2 + yk.^2)); 
kern = kern / max(kern(:));

if flag_plot
    figure()
    imagesc(kern)
    colorbar
end

end