function [pts_filled, ixs_filled] = fill_nans_with_arcs(pts, dir, flag_debug)
%FILL_NANS_WITH_ARCS Given a vector of 2D points and a 2D location for
%the center of a circle, fill with arcs the nan values of the vector. The
%arcs start and ends at the beginning and ending of the nans. The direction
%of the arc can be imposed with the dir value set to 1 (clockwise) or -1
%(anticlockwise).

if ~exist('dir','var')
    dir = 0;
end
if ~exist('flag_debug','var')
    flag_debug = false;
end

pts_filled = pts;
ixs_filled = all(isnan(pts));

if sum(ixs_filled) == 0
    error('No nans value found in the vector')
end

% Check if shift is required in case first or last elements are nans
if ixs_filled(1) > 0 || ixs_filled(end) > 0
    flag_shift = true;
else
    flag_shift = false;
end

if flag_shift
    % We re-arrange the vectors so that the first and last elements are not
    % nans. To do that we look for the first and last couple of non-nans
    % elements and shift accordingly the points.
    ixs_two_non_nan = ~ixs_filled & [~ixs_filled(2:end), ~ixs_filled(1)];
    ix_left = find(ixs_two_non_nan,1,'first');
    ix_right = find(ixs_two_non_nan,1,'last') + 1;
    if isempty(ix_left)
        error('At least two non-nan consecutive elements must be present in the vector')
    end
    pts_left = pts(:, ix_left + 1:ix_right - 1);
    pts_right = [pts(:, ix_right:end), pts(:, 1:ix_left)];
    pts_shifted = [pts_left, pts_right];
else
    pts_shifted = pts;
end

% Shift vector
ixs_shifted = all(isnan(pts_shifted));
ixs_nans_start_shifted = find(diff([0, ixs_shifted])>0);
ixs_nans_end_shifted = find(diff([ixs_shifted, 0])<0);

nf = length(ixs_nans_end_shifted);

for ix = 1:nf
    % Get indexes
    ix_start = ixs_nans_start_shifted(ix);
    ix_end = ixs_nans_end_shifted(ix);
    nh = ix_end - ix_start + 2;
    % Find limits
    R_start = norm(pts_shifted(:,ix_start - 1));
    R_end = norm(pts_shifted(:,ix_end + 1));
    th_start = atan2(pts_shifted(2,ix_start - 1), pts_shifted(1,ix_start - 1));
    th_end = atan2(pts_shifted(2,ix_end + 1), pts_shifted(1,ix_end + 1));
    if dir == 1 && th_end > th_start
        % hth must be negative
        th_start = 2*pi + th_start;
    elseif dir == -1 && th_start > th_end
        % hth must be positive
        th_end = 2*pi + th_end;
    end
    % Compute intervals
    hR = (R_end - R_start)./nh;
    hth = (th_end - th_start)./nh;
    % Sample arc
    Rvec = linspace(R_start + hR, R_end - hR, nh - 1);
    thvec = linspace(th_start + hth, th_end - hth, nh - 1);
    % Assign values
    pts_to_fill = Rvec.* [cos(thvec); sin(thvec)];
    pts_shifted(:, ix_start:ix_end) = pts_to_fill;
end

if flag_debug
    ixs_first_nan = ixs_nans_start_shifted(1);
    ixs_last_nan = ixs_nans_end_shifted(end);
    figure(), grid on, hold on, axis equal
    plot(pts_shifted(1,~ixs_shifted), pts_shifted(2,~ixs_shifted),'o')
    plot(pts_shifted(1,ixs_first_nan), pts_shifted(2,ixs_first_nan),'s-','MarkerSize',50)        
    plot(pts_shifted(1,ixs_last_nan), pts_shifted(2,ixs_last_nan),'*-','MarkerSize',50)  
    plot(pts_shifted(1,:), pts_shifted(2,:),'--')
    legend('Points','First','Last','Filling Arcs')
end

if flag_shift
    % Re-arrange the vector to the original configuration
    pts_filled(:, 1:ix_left) = pts_shifted(:, end-ix_left+1:end);
    pts_filled(:, ix_left+1:ix_right-1) = pts_shifted(:, 1:ix_right-ix_left -1);
    pts_filled(:, ix_right:end) = pts_shifted(:, ix_right-ix_left:end-ix_left);
else
    pts_filled(:, ixs_filled) = pts_shifted(:, ixs_filled);
end

end

