function [xyDist, rcoeff, pxcoeff, pycoeff] = brown_distortion(xyUndist, distortion_parameters, decentering_parameters)   
% x, y — Undistorted pixel locations. x and y are in normalized image coordinates. 
% Normalized image coordinates are calculated from pixel coordinates by translating to 
% the optical center and dividing by the focal length in pixels. Thus, x and y are dimensionless.
% k1, k2, and k3 — Radial distortion coefficients of the lens.
% p1 and p2 — Tangential distortion coefficients of the lens.

if isempty(distortion_parameters) && isempty(decentering_parameters)
    xyDist = xyUndist;
else
    k = zeros(1, 3);
    k(1:length(distortion_parameters)) = distortion_parameters;
    p = zeros(1, 2);
    p(1:length(decentering_parameters)) = decentering_parameters;
    
    rxy2 = xyUndist(1,:).^2 + xyUndist(2,:).^2;
    rcoeff = 1 + k(1)*rxy2 + k(2)*rxy2.^2 + k(3)*rxy2.^3;
    pxcoeff = 2*p(1)*xyUndist(1,:).*xyUndist(2,:) + p(2)*(rxy2 + 2*xyUndist(1,:).^2);
    pycoeff = 2*p(2)*xyUndist(1,:).*xyUndist(2,:) + p(1)*(rxy2 + 2*xyUndist(2,:).^2);
    
    xyDist = rcoeff.*xyUndist + [pxcoeff; pycoeff];
end

end