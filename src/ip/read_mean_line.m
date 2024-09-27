function data_filtered = read_mean_line(data, px_start, px_step, N_steps, M_filter, dir)

if ~exist('dir','var')
    dir = 1;
end

if dir == 1
data_filtered = nan(M_filter, size(data, 2));
for ix = 1:M_filter
    px_start_temp = px_start + ix - 1;
    data_filtered(ix,:) = mean(data(px_start_temp:px_step*2:N_steps*px_step,:), 1);
end
elseif dir == 2
data_filtered = nan(size(data, 1), M_filter);
for ix = 1:M_filter
    px_start_temp = px_start + ix - 1;
    data_filtered(:, ix) = mean(data(:, px_start_temp:px_step*2:N_steps*px_step), 2);
end 

end