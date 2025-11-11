function [ixs_occluded, Hrays, Rrays] = find_sphere_occlusions(Rsph, pos_t2p, dir_r2p, dem, ...
                                        nrs, rays_span_angle, nworkers, method_intersection, method_sampling, ...
                                        rays_step_threshold, rays_dist_first)
%FIND_SPHERE_OCCLUSIONS Find the indexes of the occluded points in a sphere
%of given radius Rsph, considering the point coordinates wrt sphere center pos_t2p, 
%the directions of the ray at each point dir_r2p and a Digital-Elevation-Map 
%in GriddedInterpolant format dem.
%The occlusions are found by sampling nrs rays from the point along the
%direction specified at each sector and computing the distance to the surface
%along the radial direction. If the distance is larger than 0 for every point
%along the ray, the sector from where the ray originated is not occluded.
%The rays are propagated up to a maximum span angle rays_span_angle.

if ~isa(dem, 'griddedInterpolant')
    error('Please provide dem as a gridded interpolant')
end

% Ray settings
rays_dist_max = Rsph(1)*tan(rays_span_angle);
if ~exist('rays_dist_first', 'var')
    rays_dist_first = 0.5*Rsph(1)*pi/min(size(dem.Values))/2;   % half the span of a pixel of DEM map
end
if ~exist('rays_step_threshold', 'var')
    rays_step_threshold = rays_dist_first/400;   % 1/400 of the first step
end

nhl = size(dir_r2p, 2);
if nhl == 0
    ixs_occluded = false;
    Hrays = [];
    Rrays = [];
    return
end

switch method_intersection

    case 'sampling'
        
        % Capping nrs to avoid memory issues
        txt = ['...', method_intersection, ' intersection, ', method_sampling, ' sampling up to ', num2str(nrs), ' ray steps in up to ', num2str(rad2deg(max(rays_span_angle))), ' deg...'];
        switch method_sampling
            case 'logspace'
                fprintf(txt)
                if numel(rays_dist_first) == 1 && numel(rays_dist_max) == 1
                    rays_dist = logspace(log10(rays_dist_first), log10(rays_dist_max), nrs)';
                else
                    rays_dist = rays_dist_first + (rays_dist_max - rays_dist_first) .* logspace(0, 1, nrs)';
                end
            case 'linspace'
                fprintf(txt)
                if numel(rays_dist_first) == 1 && numel(rays_dist_max) == 1
                    rays_dist = linspace(rays_dist_first, rays_dist_max, nrs)';
                else
                    rays_dist = rays_dist_first + (rays_dist_max - rays_dist_first) .* linspace(0, 1, nrs)';
                end
            otherwise
                error('Sampling method not recognized. use logspace or linspace.')
        end

        if nworkers > 1
            % For each ray step, find height and azimuth/elevation coordinates
            Hrays = zeros(nrs, nhl);
            Azrays = zeros(nrs, nhl);
            Elrays = zeros(nrs, nhl);
            parfor (ix = 1:nrs, nworkers)
            %for ix = 1:nrs
                pos_t2r = pos_t2p - dir_r2p.*rays_dist(ix, :);
                sph_t2r = sph_coord_fast(pos_t2r);
                Hrays(ix, :) = sph_t2r(1,:);
                Azrays(ix, :) = sph_t2r(2,:);
                Elrays(ix, :) = sph_t2r(3,:);
            end
            % Interpolate body radius at the location of the rays in IAU frame
            Rrays = find_triaxial_radius(Azrays, Elrays, Rsph) + dem(Elrays, Azrays);
            rays_step = Hrays - Rrays;
            ixs_occluded = any(rays_step <= rays_step_threshold, 1);
        else
            % Fully vectorized version
            pos_t2r = pos_t2p - dir_r2p.*reshape(rays_dist, 1, size(rays_dist, 2), nrs);
            sph_t2r = sph_coord_fast(pos_t2r);
            Hrays = permute(sph_t2r(1,:,:), [2 3 1]);
            Azrays = permute(sph_t2r(2,:,:), [2 3 1]);
            Elrays = permute(sph_t2r(3,:,:), [2 3 1]);
            Rrays = find_triaxial_radius(Azrays, Elrays, Rsph) + dem(Elrays, Azrays);
            rays_dist = Hrays - Rrays;
            ixs_occluded = any(rays_dist <= rays_step_threshold, 2)';
        end

    case 'iterative'
        
        txt = ['...', method_intersection, ' intersection, up to ', num2str(nrs), ' ray steps in up to ', num2str(rad2deg(max(rays_span_angle))), ' deg...'];
        fprintf(txt)
        if nworkers <= 1
            %---Single-call raymarching
            [ixs_occluded, iter, np_active] = raymarching(pos_t2p, dir_r2p, Rsph, dem, rays_dist_first, rays_dist_max, rays_step_threshold, nrs, false, false);
            percActiveRays = np_active/nhl;
            percRaysDone = iter/nrs;
        else
            %---Multiple-call raymarching
            % Creating pools
            idx_Pools = ceil(linspace(0, nhl, double(nworkers + uint8(1))));
            npools = length(idx_Pools) - 1;
            pos_t2p_pools = cell(1, npools);
            dir_r2p_pools = cell(1, npools);
            for idx = 1:npools
                idx_temp = idx_Pools(idx) + 1:idx_Pools(idx + 1);
                pos_t2p_pools{idx} = pos_t2p(:, idx_temp);
                dir_r2p_pools{idx} = dir_r2p(:, idx_temp);
            end
            % Parfor call
            parfor (idx = 1:npools, nworkers)
            %for idx = 1:npools
                [ixs_occluded_pool{idx}, iter{idx}, np_active{idx}] = raymarching(pos_t2p_pools{idx}, dir_r2p_pools{idx}, Rsph, dem, rays_dist_first, rays_dist_max, rays_step_threshold, nrs, false, false);
            end
            % Concatenate 
            ixs_occluded = [ixs_occluded_pool{:}];
            percActiveRays = sum([np_active{:}])./nhl;
            percRaysDone = mean([iter{:}])./nrs;
        end

        if ~isempty(percActiveRays)
            fprintf([num2str(round(1e2*percRaysDone)),'%% completed rays, ',num2str(round(1e2*percActiveRays)),'%% points left unchecked.'])
            if percActiveRays > 0.1
                fprintf('\n')
                warning('Increase occlusion_rays parameter in setting.sampling to improve rendering of shadows, if needed.')
            end
        end

    otherwise
        error('Intersection method not recognized. Use sampling or iterative.')
end

end