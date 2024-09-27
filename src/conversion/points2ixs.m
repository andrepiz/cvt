function ixs = points2ixs(pts, sz)
% Transform a vector of 2xN coordinates in (U,V) frame to linear indexes
% of the corresponding pixels of image of size sz

px_coords = round(pts + 0.5);
ixs = sub2ind(sz, px_coords(2,:), px_coords(1,:));

end