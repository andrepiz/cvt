function noise = noise_dark(ec_in, tExp, dc, sigma, seed)
% Dark Current is the collection of electrons in a pixel even when no light 
% is entering the detector. Dark current is a function of temperature, 
% typically doubling each 7 to 10 degrees temperature increment, and
% exposure time. Typical figure at room temperature (20 degrees) is 
% DC = 590 electrons/(pixel*second)
% Dark Signal Non Uniformity (DSNU) defines the different dark current 
% generation among pixels as a multiplicative factor. The effective 
% Dark Current signal on pixel i,j is obtained as DC*(1 + DSNU(i,j)).
% The DC is modelled as a Poisson distribution with a rate equal to the DC
% signal.

if exist('seed', 'var')
    if ~isempty(seed)
        rng(seed)
    end
end

if isscalar(sigma)
    dcmat = dc*(1 + sigma*randn(size(ec_in)));
elseif size(sigma, 1) == size(ec_in, 1) && size(sigma, 2) == size(ec_in, 2)
    dcmat = dc*(1 + sigma);
else
    error('Sizes of Electron Count Matrix and DSNU Standard Deviation differ')
end

dcmat = max(0, dcmat);           % no removal of electrons
noise = poissrnd_fast(dcmat).*tExp;   % dark electron noise

end

