function [ec_out, trapped, psf] = apply_blooming(ec, fwc, psf, alpha, beta, sigma_psf, tol, nMaxIter)
%APPLY_BLOOMING  Blooming model based on 2D capacity-constrained redistribution using a PSF.
%linear removal of excessive charge, there are antiblooming tension barrier
%between pixels, some pixel goes into a floating diffusion well but then
%once number increases some leak outside, other in the floating well
%   n     : HxW matrix of initial occupancy
%   M     : scalar capacity or HxW array
%   tol   : tolerance (default 1e-12)
%   maxIters : max iterations (default 1000)
%
%   Output:
%     n       : final occupancy array
%     trapped : mass that could not be placed (should be near 0)
%
% Author: Sabrina Sughi, Andrea Pizzetti — March 2026
 
if nargin < 8, nMaxIter = 1000; end  %max iteration of convergence allowed
if nargin < 7, tol = 1e-6; end %tolerance of the max value of excess of photons (when loop ocnverge)
if nargin < 6, sigma_psf = 200; end
if nargin < 5, beta = 0.02; end %offset excess over fwc where photons start the leaking
if nargin < 4, alpha = 0.05; end %percentage of excess to distribute
if nargin < 3, psf = 'gaussian'; end

% -- Matrix of photons to apply the blooming --
offset_excess = beta*fwc; %representing the floating diffusion 
n = double(ec); % Ensure doubles
% create a matrix of cut off value (fwc in this case)
if isscalar(fwc)
    fwc = fwc * ones(size(n));
else
    fwc = double(fwc);
end

% -- fwc and excess of photons to distribute with blooming --
ex = max(0, n - fwc);     % raw excess
excess = max(0, alpha*ex - offset_excess);   %excess photons to be distributed due to blooming

if isnumeric(psf)
    if isvector(psf)
        flag_1dpsf = true;
        psf = psf / sum(psf);  %to ensure conservation of spread photons
    else
        flag_1dpsf = false;
        psf = psf / sum(psf(:));  %to ensure conservation of spread photons
    end
else
    switch psf
        case 'gaussian'
            % -- Gaussian spread of excessive photons  --
            kernel_radius = ceil(3*sigma_psf); %gaussian function cut at 3*sigma
            x = -kernel_radius:kernel_radius;
            psf = exp(-x.^2/(2*sigma_psf^2));  %1D gaussian (Gaussian convolution are separable column-row)
            psf = psf / sum(psf);  %to ensure conservation of spread photons
            flag_1dpsf = true;
            kcol = psf(:);           % Nx1
            krow = psf(:).';         % 1xN

        case 'gaussian_2d'
            kernel_radius = ceil(3*sigma_psf); %gaussian function cut at 3*sigma
            x = -kernel_radius:kernel_radius;
            psf = exp(-(x.^2 + x'.^2) / (2 * sigma_psf^2));
            psf = psf / sum(psf(:));
            flag_1dpsf = false;

        otherwise
            error('PSF not supported')
    end
end

% -- Distribution of excessive photons with gaussian spread  --
%to enter the loop, photons matrix is cut at fwc value and on top of this
%the excess number of photons not eliminated with sensor architecture are
%spread 
n = n - ex; 
iter = 0;       %iteration counter initialization

%   -- Loop for excess photons convergence to tol --
while iter < nMaxIter
    iter = iter + 1;
    if flag_1dpsf
        % 2D convolution using 1D Gauss function. Gaussian 2D convolution separable row - column, 
        % this has the same result of a convolution of a 2D gaussian but way faster
        planned = conv2(conv2(excess, kcol, 'same'), krow, 'same');
    else
        % Planned redistribution using 2D convolution
        planned = conv2(excess, psf, 'same');
    end
    n = n + planned; % to previous photon matrix I add the excess photons gaussian spread
    % Compute next iteration excess
    excess = max(0, n - fwc); %remove photons higher than fwc
    % Clamp tiny numerical negatives
    n(n < 0) = 0;
    if ~any(excess(:) > tol)
        break;  % stop when the max value of excess photons is above a threshold
    end
    % Remove excess from each cell
    n = n - excess;
end
trapped = sum(max(0, n(:) - fwc(:)));  %to check the residual of excessive photons remained (should be == sum(excess(:)))

ec_out = n;

if iter == nMaxIter
    warning('Max iterations reached during application of Blooming noise. Increase nMaxIter if needed.');
end

end