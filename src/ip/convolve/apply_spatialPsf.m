function img_psf = apply_spatialPsf(img, psfs, offset_horizontal, offset_vertical, focal_length, spacing, muPixel)
%APPLY_SPATIALPSF Apply source-position-dependent, bilinearly interpolated PSFs.
%   IMG is a finite real 2-D image of integrated signal per pixel. MUPIXEL
%   and SPACING are positive image/PSF sample pitches in metres.
%   PSFS is an Nv-by-Nu cell array of finite nonnegative numeric matrices.
%   OFFSET_HORIZONTAL and OFFSET_VERTICAL give kernel-node field angles in
%   radians relative to the optical axis; positive directions are right/down.
%   FOCAL_LENGTH is the positive focal length in metres. The optical axis is
%   at the image centre; field angles follow atan(focal-plane offset / f).
%   PSFS{v,u} corresponds to these vertical/horizontal coordinates. Axes may
%   increase or decrease, but must be strictly monotonic. Singleton axes
%   specify a spatially constant PSF in that direction. Outside the node
%   grid, the nearest edge PSF is used (no extrapolation).
%
%   Each kernel is normalized independently. Kernel zero displacement is
%   sample floor(size(kernel)/2)+1, i.e. (17,17) for a 32-by-32 Zemax export.
%   Kernels must already share a common orientation and physical origin.
%   The output is double and has the input size. Zero padding loses light
%   crossing the image boundary; kernels are not renormalized at edges.
%
%   Input pixels are treated as uniform illuminated squares. PSFs are first
%   conservatively rebinned to a virtual fine grid with subdivision factor
%   ceil(MUPIXEL/SPACING). Uniform source-pixel spreading and output-pixel
%   integration are folded into an equivalent detector-sampled kernel.
%   Images are convolved at their original resolution, without oversampling.
%   This is equivalent to repeating input subpixels, applying the fine PSF,
%   and summing output subpixels, including at zero-padded image boundaries.
%   Bilinear weights are evaluated at input pixel centres before spreading,
%   so the operation conserves total signal away from image boundaries.
%   Each convolution is restricted to the nonzero weight region plus a
%   kernel-radius halo, retaining light spread outside that source region.
%   Kernel contributions are computed independently with PARFOR, then added
%   in a fixed order. Parallel workers require Parallel Computing Toolbox.
%   Reuse a pool across frames; no pool is created or deleted explicitly here.
%   Call separately for each wavelength before spectral summation.

validateattributes(img, {'numeric'}, {'real','finite','2d','nonempty'}, mfilename, 'img');
validateattributes(muPixel, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'muPixel');
validateattributes(spacing, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'spacing');
validateattributes(focal_length, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'focal_length');
assert(iscell(psfs) && ismatrix(psfs) && ~isempty(psfs), ...
    'apply_spatialPsf:InvalidGrid', 'psfs must be a nonempty 2-D cell array.');
[nv, nu] = size(psfs);
[xnodes, xorder] = validate_axis(offset_horizontal, nu, 'offset_horizontal');
[ynodes, yorder] = validate_axis(offset_vertical, nv, 'offset_vertical');
assert(all(abs(xnodes)<pi/2) && all(abs(ynodes)<pi/2), ...
    'apply_spatialPsf:InvalidAngles', 'Field angles must lie strictly between -pi/2 and pi/2.');
psfs = psfs(yorder, xorder);

[nr, nc] = size(img);
x = atan(((1:nc) - (nc+1)/2) * double(muPixel) / double(focal_length));
y = atan(((1:nr) - (nr+1)/2) * double(muPixel) / double(focal_length));
wx = axis_weights(x, xnodes);
wy = axis_weights(y, ynodes);
ratio = double(muPixel) / double(spacing);
% Avoid a spurious extra subdivision from roundoff at integer pitch ratios.
os = max(1, ceil(ratio - 16*eps(max(1, ratio))));
fine_pitch = double(muPixel) / os;
img = double(img);
img_psf = zeros(nr, nc);

% Flatten in row-major grid order to retain the previous accumulation order.
kernel_list = reshape(psfs.', [], 1);
nKernels = numel(kernel_list);
contributions = cell(nKernels, 1);
row_regions = cell(nKernels, 1);
col_regions = cell(nKernels, 1);
parfor k = 1:nKernels
    v = floor((k-1)/nu) + 1;
    u = mod(k-1,nu) + 1;
    kernel = kernel_list{k};
    validateattributes(kernel, {'numeric'}, ...
        {'real','finite','2d','nonempty','nonnegative'}, mfilename, 'PSF');
    kernel = double(kernel);
    % Scale before summing to avoid overflow for large relative values.
    peak = max(kernel(:));
    assert(peak > 0, 'apply_spatialPsf:ZeroKernel', ...
        'Each PSF must have positive total intensity.');
    kernel = kernel / peak;
    kernel = kernel / sum(kernel(:));
    active_rows = find(wy(:,v) > 0);
    active_cols = find(wx(:,u) > 0);
    if isempty(active_rows) || isempty(active_cols)
        contributions{k} = [];
        row_regions{k} = [];
        col_regions{k} = [];
        continue
    end
    kernel = rebin_kernel(kernel, double(spacing), fine_pitch);
    kernel = detector_kernel(kernel, os);
    halo = (size(kernel)-1) / 2;
    rows = max(1,active_rows(1)-halo(1)):min(nr,active_rows(end)+halo(1));
    cols = max(1,active_cols(1)-halo(2)):min(nc,active_cols(end)+halo(2));
    weights = wy(rows,v) * wx(cols,u).';
    source = img(rows,cols) .* weights;
    contributions{k} = conv2(source, kernel, 'same');
    row_regions{k} = rows;
    col_regions{k} = cols;
end
% Contributions overlap, so combine them serially without shared writes.
for k = 1:nKernels
    if isempty(contributions{k})
        continue
    end
    rows = row_regions{k};
    cols = col_regions{k};
    img_psf(rows,cols) = img_psf(rows,cols) + contributions{k};
    contributions{k} = []; % Release the cropped result after accumulation.
end
end

function [nodes, order] = validate_axis(nodes, n, name)
validateattributes(nodes, {'numeric'}, ...
    {'real','finite','vector','numel',n}, mfilename, name);
nodes = double(nodes(:));
assert(n == 1 || all(diff(nodes)>0) || all(diff(nodes)<0), ...
    'apply_spatialPsf:InvalidAxis', '%s must be strictly monotonic.', name);
[nodes, order] = sort(nodes);
end

function weights = axis_weights(points, nodes)
if numel(nodes) == 1
    weights = ones(numel(points), 1);
else
    points = min(max(points(:), nodes(1)), nodes(end));
    weights = interp1(nodes, eye(numel(nodes)), points, 'linear');
end
end

function kernel = rebin_kernel(kernel, old_pitch, new_pitch)
% Treat each original sample as integrated mass over its sample cell.
% Overlap matrices redistribute this mass without changing its total.
ry = overlap_matrix(size(kernel,1), old_pitch, new_pitch);
rx = overlap_matrix(size(kernel,2), old_pitch, new_pitch);
kernel = ry * kernel * rx.';
kernel = kernel / sum(kernel(:));
end

function kernel = detector_kernel(kernel, os)
% For a fine-grid displacement j and detector displacement d, the fraction
% of uniformly distributed source subpixels landing in detector pixel d is
% max(0, 1-abs(j-d*os)/os). Apply these triangular weights in both axes.
% Columns of each weight matrix sum to one, preserving the PSF's flux.
if os == 1
    return
end
radius = (size(kernel)-1)/2;
detector_radius = ceil(radius/os);
fine_y = -radius(1):radius(1);
fine_x = -radius(2):radius(2);
detector_y = (-detector_radius(1):detector_radius(1)).' * os;
detector_x = (-detector_radius(2):detector_radius(2)).' * os;
wy = max(0, 1-abs(detector_y-fine_y)/os);
wx = max(0, 1-abs(detector_x-fine_x)/os);
kernel = wy * kernel * wx.';
end

function overlap = overlap_matrix(n, old_pitch, new_pitch)
centres = ((1:n) - (floor(n/2)+1)) * old_pitch;
radius = max(0, ceil((max(abs(centres)) + old_pitch/2)/new_pitch - 0.5));
new_centres = (-radius:radius).' * new_pitch;
left = max(new_centres-new_pitch/2, centres-old_pitch/2);
right = min(new_centres+new_pitch/2, centres+old_pitch/2);
overlap = max(0, right-left) / old_pitch;
end
