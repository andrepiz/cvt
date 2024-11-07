function [north, east, down] = loc2ned(loc)
% Compute the north, east and down directions in a local plane tangent to a 
% perfect sphere at the specified location in a reference frame. 
% The north direction in the reference frame is assumed aligned to Z direction. 

loc = reshape(loc, 3, []);
np = size(loc, 2);

down = -loc./vecnorm(loc);
east = cross(down, repmat([0;0;1], 1, np));
east = east./vecnorm(east);
north = cross(east, down);
north = north./vecnorm(north);

end