function [coord, val, K] = project_clip_space(dir, d, projection, fov, flag_debug)

if ~exist("flag_debug","var")
    flag_debug = false;
end

posHom = [d.*dir; ones(1, length(d))];

l = min(posHom(1,:));
r = max(posHom(1,:));
b = min(posHom(2,:));
t = max(posHom(2,:));
f = max(posHom(3,:));
n = min(posHom(3,:));

if isnumeric(projection)
    K = projection;
else
    switch projection

        case 'perspective'

            % t = n*tan(fov(2)/2);
            % b = -t;
            % r = n*tan(fov(1)/2);
            % l = -r;
            % 
            alpha = (f+n)/(n-f);
            beta = -2*f*n/(n-f);
            % N = [1 0 0 0; 0 1 0 0; 0 0 alpha beta; 0 0 -1 0];
            % H = eye(4);
            % H(1, 3) = (r+l)/(r-l);
            % H(2, 3) = (t+b)/(t-b);
            % S = eye(4);
            % S(1, 1) = 2*n/(r-l);
            % S(2, 2) = 2*n/(t-b);
            % K = N*S*H;
            K = [1/tan(fov(1)/2),   0,                     0,            0;
                   0,               1/tan(fov(2)/2),       0,            0;
                   0,               0,                     alpha         beta;
                   0,               0,                    -1,            0];

        case 'orthographic'
    
            K = [2/(r-l)   0         0       -(r+l)/(r-l); ...
                 0         2/(t-b)   0       -(t+b)/(t-b);...
                 0         0         2/(f-n) -(f+n)/(f-n);...
                 0         0         0        1];

        otherwise
            error('Projection not recognized. Use perspective or orthographic or input the K you want to use.')
    end
end

p = K*posHom;
p = p./p(end,:);
coord = p([1:2], :);
val = p(3,:);

if flag_debug
    ixsPlot = round(linspace(1, length(val), 1e5));
    figure()
    grid on, hold on, axis equal
    % set(gca(), 'YDir','reverse')
    scatter(p(1,ixsPlot), p(2,ixsPlot), [], val(ixsPlot))
    colorbar
end

end