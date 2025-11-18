function [ixs_occluded, iter, np_active] = ...
    raymarching(pos_t2p, dir_r2p, Rsph, dem, rays_dist_first, ...
        rays_dist_max, rays_step_threshold, nmaxiter, flag_show_iteration, flag_debug)

nhl = size(dir_r2p, 2);
rays_dist     = rays_dist_first * ones(1, nhl);
ixs_occluded  = false(1, nhl);

if flag_debug
    figure()
    ax1 = subplot(1,2,1); grid on, hold on; title('Ray Span'); yline(rays_dist_max);
    ax2 = subplot(1,2,2); grid on, hold on; yline(rays_step_threshold); title('Ray SDF');
end

inds_act = 1:nhl;
for iter = 1:nmaxiter

    % Stop if no active rays remain
    if isempty(inds_act)
        break;
    end

    % Ray-march for each active ray
    pos_t2r = pos_t2p(:, inds_act) - dir_r2p(:, inds_act) .* rays_dist;
    sph_t2r = sph_coord_fast(pos_t2r);
    Hrays = sph_t2r(1,:);
    Azrays = sph_t2r(2,:);
    Elrays = sph_t2r(3,:);

    % Interpolate body radius at ray positions
    Rrays = find_triaxial_radius(Azrays, Elrays, Rsph) + dem(Elrays, Azrays);
    rays_sdf = Hrays - Rrays;

    % Update status arrays
    if numel(rays_dist_max) == 1
        ixs_inf = rays_dist >= rays_dist_max;
    else
        % account for different maximum distance thresholds
        ixs_inf = rays_dist >= rays_dist_max(inds_act);
    end
    ixs_occl = rays_sdf <= rays_step_threshold;
    ixs_act = ~(ixs_inf | ixs_occl);

    % Step ray distance for remaining active rays
    rays_dist = rays_dist(ixs_act) + abs(rays_sdf(ixs_act));

    % Only process currently active rays using index vector
    ixs_occluded(inds_act(ixs_occl)) = true;
    inds_act(~ixs_act) = [];

    if flag_show_iteration && mod(iter, max(1, round(nmaxiter/10))) == 0
        fprintf('%d%%: %d occluded, %d active\n', ...
            round(iter/nmaxiter*100), sum(ixs_occluded), numel(inds_act));
    end

    if flag_debug && (mod(iter, 10) == 0)
        plot(ax1, sort(rays_dist)), plot(ax2, sort(rays_sdf))
        drawnow;
    end
end

np_active = numel(inds_act);

end
