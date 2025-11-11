function rc = uv2rc(uv)
%UV2RC converts U and V coordinates to decimal row and column coordinates, 
%given U = 1 at decimal column = 0.5 and V = 1 at decimal row = 0.5.

rc = [uv(2, :) - 0.5; uv(1,:) - 0.5];

end

