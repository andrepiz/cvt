function [active, infov, observable, lit] = check_visibility(fov, N, Vinc, Vref, flag_fast)
% Check if a point observed from a camera is active, i.e. it must
% respect these three conditions:

active = false(1, size(N, 2));
[observable, lit] = deal(active);

% 1. Be inside the FOV
switch length(fov)
    case 1
        anglebearing = acos(-Vref(3,:));
        infov = anglebearing <= fov/2;
    case 2
        anglebearing = abs(atan(Vref([1 2],:)./Vref(3,:)));
        infov = all(anglebearing <= reshape(fov, 2, 1)/2);
    otherwise
        error('FOV dimension must be either 1 or 2')
end
active = infov;
if all(~active) && flag_fast
    return
end

% 2. Its surface must be visible, i.e. the angle between normal and reflection direction must be lower than 90°
delta_ref = acos(dot(N, Vref));
observable = delta_ref <= pi/2;
active = active & observable;
if all(~active) && flag_fast
    return
end

% 3. Its surface must be illuminated, i.e. the angle between normal and incident direction must be lower than 90°
delta_inc = acos(dot(N, Vinc));
lit = delta_inc <= pi/2;
active = active & lit;
if all(~active) && flag_fast
    return
end

end
