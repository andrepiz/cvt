function p = plotFov(posCam_CAM, dcm_REF2CAM, fov, h, b, col, alpha)
% h: height
% b: base

switch length(fov)
    case 1
        p = plotFovCone(posCam_CAM, dcm_REF2CAM, fov, h, b, col, alpha);
    case 2
        p = plotFovFrustum(posCam_CAM, dcm_REF2CAM, fov, h, col, alpha);
    otherwise
        error('FOV should have maximum two elements')
end

end