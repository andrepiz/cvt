function horizon = height2horizon(longrid, latgrid, height, params)
% Convert a height/displacement map to a critical horizon map. The critical
% horizon is the maximum elevation angle over which there are no 
% occlusions along any azimuthal direction for a point on a sphere.
%
% INPUTS:
% longrid [u, v]                        Grid of longitude values 
% latgrid [u, v]                        Grid of latitude values 
% height [u, v]                         2D matrix of height values at each
%                                       longitude and latitude point 
% params.granularity [1]                Granularity factor for elevation computation
%
% OUTPUTS:
% horizon [u, v]                        Map of critical horizon 


[res_v, res_u] = size(height);
lon_lims = [min(longrid(:)), max(longrid(:))];
lat_lims = [min(latgrid(:)), max(latgrid(:))];

[Finterp_height, hlon, hlat] = map2griddedInterpolant(height, lon_lims, lat_lims);
spanmax = acos(min(height(:))./max(height(:)));

horizon = pi/2*ones(res_v, res_u);
c = 0;
for ii = 1:res_u
    for jj = 1:res_v

        c = c + 1;

        lat0 = latgrid(jj, ii); 
        lon0 = longrid(jj, ii);

        lonvec_temp = lon0 - spanmax/2 : params.granularity*hlon : lon0 + spanmax/2;
        latvec_temp = lat0 - spanmax/2 : params.granularity*hlat : lat0 + spanmax/2;
        [longrid_temp, latgrid_temp] = meshgrid(lonvec_temp, latvec_temp);
        
        longrid_temp(latgrid_temp > pi/2) = pi + longrid_temp(latgrid_temp > pi/2);
        latgrid_temp(latgrid_temp > pi/2) = pi - latgrid_temp(latgrid_temp > pi/2);
        longrid_temp = wrapToPi(longrid_temp);
        
        el_vec = height2el(Finterp_height(latgrid_temp(:), longrid_temp(:)), longrid_temp(:), latgrid_temp(:), Finterp_height(lat0, lon0), lon0, lat0);
        
        horizon(jj, ii) = max(el_vec);

        disp([num2str(1e2*c/(res_u*res_v)), '%'])
    end

end

end

