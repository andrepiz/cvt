function noise = noise_prnu(ec_in, sigma, seed)
% -	Photo Response Non Uniformity
% It defines the different photoresponse among pixels. It is characterized by a Gaussian distribution. 
% Let N the mean number of electrons collected on pixel due to light effect, the PRNU is modelled
% with average value equal to 1 and a standard deviation defined by sigma. 
% The effective number of electrons on pixel i,j is obtained multiplying N*PRNU(i,j).
% Typical PRNU figure is 2% of average value (sigma_PRNU=0.02).

if exist('seed', 'var')
    if ~isempty(seed)
        rng(seed)
    end
end

if isscalar(sigma)
    matstd = sigma*randn(size(ec_in));
elseif size(sigma, 1) == size(ec_in, 1) && size(sigma, 2) == size(ec_in, 2)
    matstd = sigma;
else
    error('Sizes of Electron Count Matrix and PRNU Standard Deviation differ')
end

noise = ec_in.*matstd;

end

