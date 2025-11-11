function cart_coord = cart_coord_fast(sph_coord)
%convert spherical coordinate in cartesian coordinates
if ~(size(sph_coord,1)==3 )
    error('Input must be 3 x N')
end
if ~(size(sph_coord,3)==1 )
    error('Use cart_coord for arrays with a 3rd dimension')
end

r = sph_coord(1,:);
az = sph_coord(2,:);
el = sph_coord(3,:);

rcosel = r.*cos(el);
rsinel = r.*sin(el);

cart_coord  = [ rcosel.*cos(az)
                rcosel.*sin(az)
                rsinel];

end
