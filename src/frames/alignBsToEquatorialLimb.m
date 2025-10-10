function [q_CAMI2CAM, uLimb, vLimb, pitch, roll, bodyAngSize, bodyPxSize] = alignBsToEquatorialLimb(R, d, phase_angle, f, muPixel, res_px, bsAngle)
% Starting from the boresight yaw angle (angle around boresight), align the
% pitch and yaw axis so that the intersection of the body limb with its
% equator is at the optical center.

bodyAngSize = 2*asin(R./d);
bodyPxSize = 2*tan(bodyAngSize/2).*f./muPixel;

yaw = bsAngle;
uLimb = res_px(1)/2 + bodyPxSize/2.*cos(yaw);
vLimb = res_px(2)/2 - bodyPxSize/2.*sin(yaw);

pitch = -sign(phase_angle).*atan((uLimb-res_px(1)/2).*muPixel./f);
roll = sign(phase_angle).*atan((vLimb-res_px(2)/2).*muPixel./f);

q_CAMI2CAM = euler_to_quat([yaw; pitch; roll], struct('sequence','zyx'));

end