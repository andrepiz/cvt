function per = get_fov_perimeter(fov, nper, flag_debug)
%GET_FOV_PERIMETER Find the FOV perimeter as the space of LOS sampled with
%specified number of points nper.
%The FOV is the full-span angle of the horizontal (first element) and 
%vertical (second element) directions. If a single value is provided, 
% the FOV is a cone, if two values are provided, the FOV is a frustrum.
%Boresight of camera is assumed aligned to +Z direction
bs = [0; 0; 1];

if ~exist('flag_debug','var')
    flag_debug = false;
end

zvec = zeros(1, nper);
switch length(fov)
    case 1 
        lvec1 = linspace(1, 0, nper);
        ovec = ones(1, nper);
        roll_vec = 0.5*fov*ovec;
        pitch_vec = zvec;
        yaw_vec = 2*pi*lvec1;
        q_bs2fov = euler_to_quat([roll_vec; pitch_vec; yaw_vec]);
        per = rotframe(bs, q_bs2fov, true);

    case 2
        nper = 8*(ceil(nper/8));
        ovec = ones(1, nper/4);
        lvec1_left = -flip(linspace(8/nper, 1, nper/8).^cos(fov(1)/2));
        lvec1_right = linspace(0, 1 - 8/nper, nper/8).^cos(fov(1)/2);
        lvec1 = [lvec1_left, lvec1_right];
        lvec2_left = -flip(linspace(8/nper, 1, nper/8).^cos(fov(2)/2));
        lvec2_right = linspace(0, 1 - 8/nper, nper/8).^cos(fov(2)/2);
        lvec2 = [lvec2_left, lvec2_right];
        fun_sampling = @(x) sign(x).*sqrt(sin(x).^2./(1-sin(x).^2));
        x = [fun_sampling(0.5*fov(1)*lvec1),  fun_sampling(0.5*fov(1)*ovec),  -fun_sampling(0.5*fov(1)*lvec1),  -fun_sampling(0.5*fov(1)*ovec)];
        y = [fun_sampling(0.5*fov(2)*ovec),  fun_sampling(-0.5*fov(2)*lvec2),  -fun_sampling(0.5*fov(2)*ovec),  fun_sampling(0.5*fov(2)*lvec2)];
        z = ones(1, nper);
        per = [x; y; z];
        per = per./vecnorm(per);

    otherwise
        error('FOV should have either 1 or 2 dimensions')
end

       
if flag_debug
    zvec = zeros(1, nper);
    figure(), grid on, hold on, axis equal
    quiver3(zvec, zvec, zvec, per(1, :), per(2, :), per(3, :))
    xlabel('x'), ylabel('y'), zlabel('z'), view([10, 20])
    title('FOV')
end

end

