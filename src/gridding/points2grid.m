function [Zidw, Zmin, Zmax, densitymap] = points2grid(pointcloud, weighter, dx, flag_plot)

dy = dx;
%search_radius = 20;
search_radius = sqrt(2).*dx;

x = pointcloud(:,1)-min(pointcloud(:,1));
y = pointcloud(:,2)-min(pointcloud(:,2));
z = pointcloud(:,3);

minx = min(x);
maxx = max(x);
miny = min(y);
maxy = max(y);
minz = min(z);
maxz = max(z);
xx = (minx+dx):dx:(maxx);
yy = (miny+dy):dy:(maxy);

[X,Y] = meshgrid(xx,yy');

for j=1:length(yy)
for k=1:length(xx)
%plot(X(j,k), Y(j,k), 'bo')
%Use this to search only on those points within the actual grid cell (within dx and dy of the node)
%tf = x<= X(j,k)+dx./2 & x>=X(j,k)-dx./2 & y<=Y(j,k)+dy./2 & y>=Y(j,k)-dy./2;
%locs=find(tf);
%if length(locs)==0
% Zmin(j,k)=NaN;
% Zmean(j,k)=NaN;
% Zmax(j,k)=NaN;
%else
% Zmin(j,k)=min(z(locs));
% Zmean(j,k)=mean(z(locs));
% Zmax(j,k)=max(z(locs));
%end
tf = x<= X(j,k)+search_radius & x>=X(j,k)-search_radius & y<=Y(j,k)+search_radius & y>=Y(j,k)-search_radius;
locs=find(tf);
localx = x(locs);
localy = y(locs);
localz = z(locs);
dist = sqrt((localx-X(j,k)).^2 + (localy-Y(j,k)).^2);
locs_radius=find(dist<=search_radius);
plot(localx(locs_radius), localy(locs_radius), 'm.')
axis([minx maxx miny maxy])
title('points and nodes')
if isempty(locs_radius)
Zmin(j,k)=NaN;
Zmean(j,k)=NaN;
Zmax(j,k)=NaN;
densitymap(j,k)=NaN;
else
Zmin(j,k)=min(localz(locs_radius));
Zmean(j,k)=mean(localz(locs_radius));
Zmax(j,k)=max(localz(locs_radius));
Zidw(j,k) = sum(localz./(dist.^weighter))./sum(1./(dist.^weighter));
densitymap(j,k)=length(localz);
end
end
end

if flag_plot
figure(1)
clf

subplot(2,3,1)
plot(x,y, 'k.')
hold on
plot(X, Y, 'r+')
xlabel('Easting (ft)')
ylabel('Northing (ft)')

clims = [0 260];
subplot(2,3,2)
plot(minx, miny, 'k.')
hold on
imagesc(xx', yy, Zmax, clims)
axis([minx maxx miny maxy minz maxz])
%colorbar
title('Zmax')
subplot(2,3,3)
plot(minx, miny, 'k.')
hold on
imagesc(xx', yy, Zmin, clims)
axis([minx maxx miny maxy minz maxz])
%colorbar
title('Zmin')
subplot(2,3,4)
plot(minx, miny, 'k.')
hold on
imagesc(xx', yy, Zmean, clims)
axis([minx maxx miny maxy minz maxz])
%colorbar
title('Zmean')
subplot(2,3,5)
plot(minx, miny, 'k.')
hold on
imagesc(xx', yy, Zidw, clims)
axis([minx maxx miny maxy minz maxz])
%colorbar
title('Zidw')
subplot(2,3,6)
plot(minx, miny, 'k.')
hold on
imagesc(xx', yy, Zmax-Zmin, clims)
axis([minx maxx miny maxy minz maxz])
colorbar
title('Zdif')
figure(2)
clf
surfl(X,Y,Zmax)
shading interp
hold on
surfl(X,Y,Zmin)
colormap(jet)
xlabel('Easting (ft)')
ylabel('Northing (ft)')
zlabel('Elevation (ft)')

figure(3)
clf
plot(minx, miny, 'k.')
hold on
imagesc(xx', yy, densitymap)
axis([minx maxx miny maxy])
colormap cool
colorbar
xlabel('Easting (ft)')
ylabel('Northing (ft)')
title('point density per pixel')
end
end