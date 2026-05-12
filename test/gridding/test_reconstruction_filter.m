%% Test script for image reconstruction filters
cvt_install()

% Load data
%load('test_gridding_sphere1M.mat')
load('test_gridding_moonlimb3M.mat')

flag_parallelization = true;
flag_scaling = true;
nworkers = 14;
% histweight algorithm
gridding_algorithm = 'area';
% reconstruction algorithm
reconstruction_granularity = 2;
kern_antialiasing = gaussianKernel(reconstruction_granularity, sqrt(reconstruction_granularity), false);
reconstruction_filter = 'bilinear';
% shiftedquantization algorithm
gridding_shift = 2;
weighting_variance = 1;
kern_weighting = gaussianKernel(gridding_shift, sqrt(weighting_variance), false);

% number of tests to do
Ntest = 8;

if flag_scaling 
    scaling_factor = max(PL_point(:))*4;
else
    scaling_factor = 4;
end

% pre-processing
vals = PL_point/scaling_factor;
img = zeros(res_px(1), res_px(2), Ntest);
c = 1;

%% HISTWEIGHT AREA
tic
if flag_parallelization
    temp = parhistweight(coords, vals, [1, res_px(1); 1 res_px(2)], 1, gridding_algorithm, nworkers);
else
    temp = histweight(coords, vals, [1, res_px(1); 1 res_px(2)], 1, gridding_algorithm, false);
end
cputime(c) = toc;
img(:,:,c) = temp;
c = c + 1;

%% UNRESAMPLED SUM
tic
temp = quantization(coords, vals, [1, res_px(1); 1 res_px(2)], 1, 'sum');
cputime(c) = toc;
img(:,:,c) = temp;
c = c + 1;

%% UPSAMPLED SUM + RECONSTRUCTION
tic
temp = quantization(coords, vals, [1, res_px(1); 1 res_px(2)], reconstruction_granularity, 'sum');
temp = downsamplingreconstruction(temp, reconstruction_granularity, [], reconstruction_filter, reconstruction_granularity);
cputime(c) = toc;
img(:,:,c) = temp;
c = c + 1;

%% UPSAMPLED SUM + ANTI-ALIASING + RECONSTRUCTION
tic
temp = quantization(coords, vals, [1, res_px(1); 1 res_px(2)], reconstruction_granularity, 'sum');
temp = downsamplingreconstruction(temp, reconstruction_granularity, [], reconstruction_filter, reconstruction_granularity);
cputime(c) = toc;
img(:,:,c) = temp;
c = c + 1;

%% SHIFTED UNRESAMPLED SUM
tic
temp = shiftedquantization(coords, vals, [1, res_px(1); 1 res_px(2)], 1, 'sum', gridding_shift);
cputime(c) = toc;
img(:,:,c) = temp;
c = c + 1;


%% SHIFTED WEIGHTED UNRESAMPLED SUM
tic
temp = shiftedquantization(coords, vals, [1, res_px(1); 1 res_px(2)], 1, 'sum', gridding_shift, kern_weighting);
cputime(c) = toc;
img(:,:,c) = temp;
c = c + 1;

%% SHIFTED UPSAMPLED SUM
tic
temp = shiftedquantization(coords, vals, [1, res_px(1); 1 res_px(2)], reconstruction_granularity, 'sum', gridding_shift);
temp = downsamplingreconstruction(temp, reconstruction_granularity, [], reconstruction_filter, reconstruction_granularity);
cputime(c) = toc;
img(:,:,c) = temp;
c = c + 1;

%% SHIFTED WEIGHTED UPSAMPLED SUM
tic
temp = shiftedquantization(coords, vals, [1, res_px(1); 1 res_px(2)], reconstruction_granularity, 'sum', gridding_shift, kern_weighting);
temp = downsamplingreconstruction(temp, reconstruction_granularity, [], reconstruction_filter, reconstruction_granularity);
cputime(c) = toc;
img(:,:,c) = temp;
c = c + 1;

%% VISUAL COMPARISON

timg = imtile(img, 'GridSize',[2 Ntest/2]);

figure()
imshow(timg)

imgref = repmat(img(:,:,1), 1, 1, Ntest);
imgerr = imabsdiff(imgref, img);

figure()
subplot(1,2,1)
plot(cputime)
ylabel('CPU Time [s]')
xlabel('method')
subplot(1,2,2)
plot(squeeze(sum(imgerr, [1 2])))
ylabel('Error [DN]')
xlabel('method')


%%

% fh1 = figure();
% grid on, hold on
% grid minor
% scatter(coords(2,:), coords(1,:), [], PL_point)
% colormap('parula');
% axis equal
% fh1.CurrentAxes.YDir = 'reverse';
% colormap('parula');
% col = colorbar;
% col.Label.String = 'PL [m^2]';
% xlabel('u [px]')
% ylabel('v [px]')
% xlim([0, res_px(1)])
% ylim([0, res_px(2)])