function fov = resfpx2fov(res, fpx)
%FPX2FOV return fov from resolution and pixel pitch

fov = 2*atan(res/2/fpx);

end

