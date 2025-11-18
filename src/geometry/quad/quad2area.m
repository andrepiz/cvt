function area = quad2area(lon1, lat1, lon2, lat2, dem, radius, dcm_CSF2IAU)

% Extract longitude/latitude lower limits for each sector
[lon1Grid, lat1Grid] = meshgrid(lon1, lat1);

% Extract longitude/latitude upper limits for each sector
[lon2Grid, lat2Grid] = meshgrid(lon2, lat2);

lonGrid = [lon1Grid(:), lon2Grid(:), lon2Grid(:), lon1Grid(:)];
latGrid = [lat1Grid(:), lat1Grid(:), lat2Grid(:), lat2Grid(:)];

lonGrid = lonGrid(:)';
latGrid = latGrid(:)';
RBody = find_triaxial_radius(lonGrid, latGrid, radius);

if ~isempty(dem)
    PNorm_CSF = cart_coord_fast([RBody; lonGrid; latGrid]);
    PNorm_IAU = dcm_CSF2IAU*PNorm_CSF;
    sph_IAU = sph_coord_fast(PNorm_IAU);
    RGrid = RBody + dem(sph_IAU(3,:), sph_IAU(2,:));
else
    RGrid = Rbody;
end
P = cart_coord_fast([RGrid; lonGrid; latGrid]);
np = size(P, 2);
area = vertex2area(P(:,1:np/4), P(:,np/4+1:np/2), P(:,np/2+1:3*np/4), P(:,3*np/4+1:np));


flag_debug = false;
if flag_debug

%%
figure;
hold on;
for idx = 1:1:1000
    idxs = [idx; np/4+1+idx; np/2+1+idx; 3*np/4+1+idx];
    fill3(P(1,idxs), P(2,idxs), P(3,idxs), area(idx));
end
axis equal

%%
figure;
hold on;
for idx = 1:100:108000
    plot3(squeeze(P(1,idx,1)), squeeze(P(2,idx,1)), squeeze(P(3,idx,1)), 'o');
end
axis equal
end

end