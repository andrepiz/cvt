function [ixs_occluded, Hray, Rray] = find_sphere_occlusions(Rsph, pos_t2p, dir_r2p, dem, ...
                                        nrs, rays_span_angle, nworkers, method_intersection, method_sampling)
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

switch method_intersection
    case 'sampling'
        nhl = size(dir_r2p, 2);
        % Sample the ray
        ray_span_max = Rsph(1)*tan(rays_span_angle);
        ray_span_min = Rsph(1)*pi/2/min(size(dem.Values));   % half the span of a pixel of DEM map
        switch method_sampling
            case 'logspace'
                ray_span_vec = logspace(log10(ray_span_min), log10(ray_span_max), nrs);
            case 'linspace'
                ray_span_vec = linspace(ray_span_min, ray_span_max, nrs);
        end
        Hray = zeros(nrs, nhl);
        Azray = zeros(nrs, nhl);
        Elray = zeros(nrs, nhl);
        % For each ray step, find height and azimuth/elevation coordinates
        parfor (ix = 1:nrs, nworkers)
            pos_t2r = pos_t2p - dir_r2p.*ray_span_vec(ix);
            sph_t2r = sph_coord(pos_t2r);
            Hray(ix, :) = sph_t2r(1,:);
            Azray(ix, :) = sph_t2r(2,:);
            Elray(ix, :) = sph_t2r(3,:);
        end
        % Interpolate body radius at the sector location of the ray in IAU frame
        Rray = find_triaxial_radius(Azray, Elray, Rsph) + dem(Elray, Azray);
        d_rays = Hray - Rray;
        ixs_occluded = any(d_rays <= 0, 1);
    case 'iterative'
        error('method not supported')

end