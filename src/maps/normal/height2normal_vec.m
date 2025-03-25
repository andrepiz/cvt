function [normal, debug] = height2normal_vec(longrid, latgrid, height, params)
% Convert a height/displacement map to a normal map using either a plane
% fitting algorithm or a mean normal algorithm
%
% INPUTS:
% longrid [u, v]                        Grid of longitude values 
% latgrid [u, v]                        Grid of latitude values 
% height [u, v]                         2D matrix of height values at each
%                                       longitude and latitude point 
% params.flag_debug [true/false]        Plot points and normals
% params.fitting_error_threshold [1]    If plane fitting error is larger than 
%                                       this threshold, the normal is computed as the mean of the normals of the
%                                       neighbouring triangles
% params.frame [body/local]             If body frame is selected, the normal is stored defined in the body-fixed frame
%                                       If local frame is selected, angles of the normal are computed with 
%                                       respect to the east direction (red channel), north direction (green channel) 
%                                       and up direction (blue channel) in the local frame
%                                       
% OUTPUTS:
% normal [u, v, 3]                      Map of 3d normals at each longitude and
%                                       latitude
% debug

if ~exist('params','var')
    params = struct();
end
if ~isfield(params,'flag_debug')
    params.flag_debug = false;
end
if ~isfield(params,'fitting_error_threshold')
    md = mean(abs([reshape(diff(height, [], 1), 1, []), reshape(diff(height, [], 2), 1, [])]));
    params.fitting_error_threshold = 1e-5*md;
end
if ~isfield(params,'frame')
    params.frame = 'local'; 
end

fitting_error_threshold = params.fitting_error_threshold;
frame = params.frame;
flag_debug = params.flag_debug;
if flag_debug
    fh = figure();
    ax = axes(fh);
    grid on, hold on, axis equal
    view([1 1 1])
end

% Sizes
[nlat, nlon] = size(height);
np = nlat*nlon;

% Kernel shifts
kernSize = 1;
kernComb = permn(-kernSize:kernSize, 2);
nc = size(kernComb, 1);

% Compute kernels spherical coordinates
azkern = cast(zeros(np, nc), class(height));
elkern = cast(zeros(np, nc), class(height));
rkern = cast(zeros(np, nc), class(height));
for ix = 1:nc
    % [X, Y] = [latShift, lonShift]
    % [-1, -1] means shift latitudes northward and longitude westwards.
    % 9: north-west
    % 8: north-center
    % 7: north-east
    % 6: center-west
    % 5: center_center
    % 4: center_east
    % 3: south-west
    % 2: south-center
    % 1: south-east
    height_shifted = circshift(height, kernComb(ix,:));
    longrid_shifted = circshift(longrid, [0, kernComb(ix,2)]);
    latgrid_shifted = circshift(latgrid, [kernComb(ix,1), 0]);
    azkern(:, ix) = longrid_shifted(:);
    elkern(:, ix) = latgrid_shifted(:);
    rkern(:, ix) = height_shifted(:);
end

% Compute normals
Nx = cast(zeros(1, np), class(height));
Ny = Nx;
Nz = Nx;
if flag_debug
    algo = nan(1, np);
    err = nan(1, np);
    cnd = nan(1, np);
else
    algo = [];
    err = [];
    cnd = [];
end

parfor ix = 1:np
%for ix = 1:np

    rkerntemp = rkern(ix, :);
    azkerntemp = azkern(ix, :);
    elkerntemp = elkern(ix, :);

    % Check if points of the kernel are in the same emisphere
    ixs_same_emisphere = abs(elkerntemp(5) - elkerntemp) <= 0.99*pi;
    nctemp = sum(ixs_same_emisphere);

    % Extract kernel and mid points
    cartKern = cart_coord([rkerntemp(ixs_same_emisphere); azkerntemp(ixs_same_emisphere); elkerntemp(ixs_same_emisphere)]);
    cartMid = cart_coord([rkerntemp(5); azkerntemp(5); elkerntemp(5)]);

    x = cartKern(1,:);
    y = cartKern(2,:);
    z = cartKern(3,:);

    % Construct the matrix A and vector b for Ax = b
    A = [x', y', ones(nctemp, 1)];
    b = z';
    
    % Solve the least squares problem to find the plane coefficients
    coeffs = A \ b;

    err_temp = mean(abs(A*coeffs-b))./norm(cartMid);

    if err_temp < fitting_error_threshold
        % The plane equation is: z = Ax + By + C
        % Rearrange to Ax + By - z + C = 0
        % A vector normal to the plane is [A, B, -1]
        vec_norm = [coeffs(1); coeffs(2); -1];
        algo_temp = 1;
    else
        % Compute the normal as the mean of the normals of each
        % triangle
        pts = [x; y; z];
        vec = pts - repmat(cartMid, 1, nctemp);
        vec_norm = vec2norm(vec, nctemp);
        algo_temp = 2;
    end

    % Final normalization
    vec_normalized = vecnormalize(vec_norm);

    % Invert to make it outward
    if dot(vec_normalized, vecnormalize(cartMid)) < 0
        vec_normalized = -vec_normalized;
    end

    switch frame
        case 'body'
            Nx(ix) = vec_normalized(1);
            Ny(ix) = vec_normalized(2);
            Nz(ix) = vec_normalized(3);
        case 'local'
            [north, east, down] = loc2ned(cartMid);
            Nx(ix) = dot(vec_normalized, east);
            Ny(ix) = dot(vec_normalized, north);
            Nz(ix) = dot(vec_normalized, -down);
        otherwise
            error('Frame not recognized. Use body or local options')
    end
    
    if flag_debug
        algo(ix) = algo_temp;
        err(ix) = err_temp;
        cnd(ix) = cond(A);

        % c = 1;
        % kvec = vecnorm(cartMid);
        % try
        % q = vec./vecnorm(vec);
        % for k = 1:size(A,1)
        %     pl(c) = quiver3(ax, pts(1,k), pts(2,k), pts(3,k), kvec*q(1,k), kvec*q(2,k), kvec*q(3,k), 1, 'Color','k'); 
        %     pl(c + 1)= text(ax, pts(1,k), pts(2,k), pts(3,k), num2str(k)); 
        %     c = c + 2; 
        % end
        % catch
        % end
        % camtarget(cartMid);
        % plot3(ax, cartMid(1), cartMid(2), cartMid(3), 'ro');
        % pl(c) = quiver3(ax, cartMid(1), cartMid(2), cartMid(3), kvec*vec_normalized(1), kvec*vec_normalized(2), kvec*vec_normalized(3), 2, 'Color','g');
        % drawnow
        % delete(pl)
    end
    %disp(['Progress: ',num2str(1e2*ix/np),'%'])
end

normal = cat(3, reshape(Nx, nlat, nlon), reshape(Ny, nlat, nlon), reshape(Nz, nlat, nlon));

debug.algo = algo;
debug.err = err;
debug.cond = cnd;


end

function vec_norm = vec2norm(vec, sz)
    vecn = vecnormalize(vec);
    switch sz
        case 9
            vec_cross = cross(vecn(:, [3 6 9 8 7 4 1 2]), vecn(:, [6 9 8 7 4 1 2 3]));
        case 6
            if isnan(sum(vecn(:, 2)))
                vec_cross = cross(vecn(:, [3 6 5 4]), vecn(:, [6 5 4 1]));
            elseif isnan(sum(vecn(:, 3)))
                vec_cross = cross(vecn(:, [1 2 4 6]), vecn(:, [2 4 6 5]));
            elseif isnan(sum(vecn(:, 4)))
                vec_cross = cross(vecn(:, [2 1 3 5]), vecn(:, [1 3 5 6]));
            elseif isnan(sum(vecn(:, 5)))
                vec_cross = cross(vecn(:, [4 1 2 3]), vecn(:, [1 2 3 6]));
            else
                error(['Singularity detected at point (', num2str(ii), ',',num2str(jj),')'])
            end
        case 4
            if isnan(sum(vecn(:, 1)))
                vec_cross = cross(vecn(:, [2 4]), vecn(:, [4 3]));
            elseif isnan(sum(vecn(:, 2)))
                vec_cross = cross(vecn(:, [4 3]), vecn(:, [3 1]));
            elseif isnan(sum(vecn(:, 3)))
                vec_cross = cross(vecn(:, [1 2]), vecn(:, [2 4]));
            elseif isnan(sum(vecn(:, 4)))                  
                vec_cross = cross(vecn(:, [3 1]), vecn(:, [1 2]));
            else
                error(['Singularity detected at point (', num2str(ii), ',',num2str(jj),')'])
            end
        otherwise
            error('Kernel size not allowed')
    end
    norms = vecnormalize(vec_cross);
    vec_norm = mean(norms, 2, 'omitnan');
end