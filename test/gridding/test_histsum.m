granularity = 1;
flag_from_sim = true;
winsize = 1;

%%
if flag_from_sim
    K = rend.camera.K;
    res_px = rend.camera.res_px;
    coordsXYZ = rend.cloud.coords(:, rend.cloud.ixsValid);
    dValues = rend.cloud.values(:, rend.cloud.ixsValid);
    
    coordsUV = K*coordsXYZ;
    coordsRCP = coordsUV(1:2,:)./coordsUV(3,:);
    coords = [coordsRCP(2,:); coordsRCP(1,:)]; % Coordinates must be provided in indexing frame (rows, columns, pages)
    limits = [1, res_px(1); 1 res_px(2)];
    %limits = [1, 800; 1 1024];
    values = dValues*granularity;
else
    n = 4;
    res_px = [4 4];
    coords = rand(2, n).*[res_px(1); res_px(2)];
    values = ones(1, n);
    limits = [1, res_px(1); 1 res_px(2)];
end

%%

%--SUM
tic
bins_histsum_legacy = quantization(coords, values, limits, granularity, 'method','sum');
time_histsum_legacy = toc

tic
bins_histsum = histsum_2d(coords, values, limits, granularity);
time_histsum = toc

tic
bins_parhistum = parhistsum_2d(coords, values, limits, granularity);
time_parhistsum = toc


%--WEIGHTEDSUM
% tic
% bins_histweight_legacy = histweight_2d_legacy(coords, values, limits, granularity);
% time_histweight_legacy = toc
% 
% tic
% bins_parhistweight_legacy = parhistweight_2d_legacy(coords, values, limits, granularity);
% time_parhistweight_legacy = toc

tic
bins_histweight = histweight_2d(coords, values, limits, granularity, 2, winsize);
time_histweight = toc

tic
bins_parhistweight = parhistweight_2d(coords, values, limits, granularity, 2, winsize, 1/3, false, 6);
time_parhistweight = toc

% PLOT
 
%SUM
figure(),
subplot(2,2,1), axis equal, imagesc(bins_histsum_legacy), title(['Histsum, legacy, ', num2str(time_histsum_legacy),' s'])
subplot(2,2,2), axis equal, imagesc(bins_histsum), title(['Histsum, new, ', num2str(time_histsum),' s'])
subplot(2,2,3), axis equal, imagesc(bins_parhistum), title(['Parhistsum, new, ', num2str(time_parhistsum),' s'])

%WEIGTHEDSUM
figure(), axis equal
%subplot(2,2,1), axis equal, imagesc(bins_histweight_legacy), title(['Histweight, legacy, ', num2str(time_histweight_legacy),' s'])
subplot(2,2,2), axis equal, imagesc(bins_histweight), title(['Histweight, new, ', num2str(time_histweight),' s'])
%subplot(2,2,3), axis equal, imagesc(bins_parhistweight_legacy), title(['Parhistweight, legacy, ', num2str(time_parhistweight_legacy),' s'])
subplot(2,2,4), axis equal, imagesc(bins_parhistweight), title(['Parhistweight, new, ', num2str(time_parhistweight),' s'])

%% ERRORS
disp(['Relative error histsum: ', num2str(1e2*sum(imabsdiff(bins_histsum_legacy, bins_histsum)./sum(bins_histsum_legacy(:)),'all')),'%'])
disp(['Relative error parhistsum: ', num2str(1e2*sum(imabsdiff(bins_histsum_legacy, bins_parhistum)./sum(bins_histsum_legacy(:)),'all')),'%'])
disp(['Relative error histweight: ', num2str(1e2*sum(imabsdiff(bins_histweight_legacy, bins_histweight)./sum(bins_histweight_legacy(:)),'all')),'%'])
fprintf(['\nRelative error parhistweight: ', num2str(1e2*sum(imabsdiff(bins_parhistweight_legacy, bins_parhistweight)./sum(bins_histweight_legacy(:)),'all')),'%'])