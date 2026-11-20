function V = rotateGriddedInterpolation(F, lonA, latA, dcm_A2B)
%rotateGriddedInterpolation Interpolate an B-frame field at points given in A frame
% using a griddedInterpolant F defined in A frame.
%
%   V = rotateGriddedInterpolation(F, lonA, latA, dcm_A2B)
%
%   F is a griddedInterpolant defined over a longitude/latitude grid in the
%   B (body-fixed) frame: F.GridVectors = {lonB, latB} (degrees),
%   ndgrid-style. lonA/latA are query points (degrees, any size, same
%   size as each other) expressed in the A frame. This function rotates
%   each query direction from A into B using dcm_A2B, then samples F
%   there, returning the value(s) at those physical directions.
%
%   INPUTS
%       F         - griddedInterpolant with F.GridVectors = {lonB, latB}
%                   and F.Values sized [numel(lonB) x numel(latB)].
%       lonA    - longitude(s) of query point(s) in A frame (deg), any
%                   array size.
%       latA    - latitude(s) of query point(s) in A frame (deg), same
%                   size as lonA.
%       dcm_A2B - rotation taking vectors expressed in A to B
%
%   OUTPUT
%       V         - interpolated value(s) of F at the given A points,
%                   same size as lonA/latA.

sz = size(lonA);
N  = numel(lonA);

% --- Rotate from A into B, working with flattened row vectors --------
dir_A = cart_coord_fast([ones(1,N); lonA(:).'; latA(:).']);   % 3 x N
dir_B = dcm_A2B * dir_A;                                       % 3 x N
sph_B = sph_coord_fast(dir_B);                                 % 3 x N -> [r; lat; lon]

% --- Query the B-frame interpolant, then reshape back to input shape -
V = reshape(F(sph_B(3,:), sph_B(2,:)), sz);

end