function noise = noise_shot(ec_in, seed)
% Shot noise changes the number of photons hitting a pixel during the exposure
% with a statistical fluctaution. According to the laws of quantum mechanics, 
% the probability is Poisson distributed with rate factor equal to the
% collected number of photons. We use electrons instead of photons as 
% there is only a scaling term from photons which is the QE.

if exist('seed', 'var')
    if ~isempty(seed)
        rng(seed)
    end
end

ec_out = poissrnd(ec_in);

noise = ec_out - ec_in;

end

