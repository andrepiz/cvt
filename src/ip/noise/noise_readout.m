function noise = noise_readout(ec_in, fwc, sigma, seed)
% -	Readout noise (temporal noise)
% Readout noise is assumed Gaussian distribution, mean value equal to zero.
% Typical figure of standard deviation is 100 electrons for a 74000 FWC,
% about 2%.
% EMVA1288 Standard reports about 2.5 electrons CMOS and 8-10 electrons
% CCD

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
    error('Sizes of Electron Count Matrix and Readout Standard Deviation differ')
end

noise = matstd;
% can't read less than 0 or more than FWC
noise(ec_in + noise > fwc) = 0;
noise(ec_in + noise < 0) = 0;

end

