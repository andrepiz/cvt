function per = get_fov_perimeter(fov, nper, flag_debug)
%GET_FOV_PERIMETER Find the FOV perimeter as the space of LOS sampled with
%specified number of points nper.
%Boresight of camera is assumed aligned to +Z direction
bs = [0; 0; 1];

if ~exist('flag_debug','var')
    flag_debug = false;
end

zvec = zeros(1, nper);
switch length(fov)
    case 1
        lvec = linspace(0, 1, nper);
        ovec = ones(1, nper);
        roll_vec = fov*ovec;
        pitch_vec = zvec;
        yaw_vec = 2*pi*lvec;
        q_bs2fov = euler_to_quat([roll_vec; pitch_vec; yaw_vec]);
        per = rotframe(bs, q_bs2fov, true);

    case 2
        nper = 4*(ceil(nper/4));
        ovec = ones(1, nper/4);
        lvec = linspace(-1, 1 - 4/nper, nper/4);
        x = [sin(fov(1)*lvec), ovec*sin(fov(1)), -sin(fov(1)*lvec), -ovec*sin(fov(1))];
        y = [-ovec*sin(fov(2)), sin(fov(2)*lvec), ovec*sin(fov(2)), -sin(fov(2)*lvec)];
        z = ones(1, nper);
        per = [x; y; z];
        per = per./vecnorm(per);

    otherwise
        error('FOV should have either 1 or 2 dimensions')
end

       
if flag_debug
    figure(), grid on, hold on, axis equal
    quiver3(zvec, zvec, zvec, per(1, :), per(2, :), per(3, :))
    xlabel('x'), ylabel('y'), zlabel('z'), view([10, 20])
    title('FOV')
end

end

