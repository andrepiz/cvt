function fh = plot_frames_and_vectors(R_frames2ref, R_pos_ref, v_ref, v_pos_ref, fh, R_lg, v_lg, R_sz, R_lw, v_lw, col)
    
if ~exist('fh','var')
    fh = figure();
    grid on
    hold on
    axis equal
    view([1 1 1])
end

[a,b,n] = size(R_frames2ref);
[c,m] = size(v_ref);

if ~exist('R_lg','var')
    R_lg = cellstr(num2str([1:n]'));
end

if ~exist('v_lg','var')
    v_lg = cellstr(num2str([1:m]'));
end

if ~exist('R_sz','var')
    R_sz = ones(1, n);
end

if ~exist('R_lw','var')
    R_lw = ones(1, n);
end

if ~exist('v_lw','var')
    v_lw = ones(1, m);
end

if ~exist('col','var')
    % random colors
    %col = rand(m+n, 3);
    col = get(gca, 'ColorOrder');
end

if a~=3 || b~=3 && n~=0
    error('Provide R as [3x3xN]')
end

if c~=3 && m~=0
    error('Provide v as [3xM]')
end

col = interp1(linspace(0, 1, size(col,1)), col, linspace(0, 1, n+m+1));
for i = 1:n
    c = col(i,:);
    R = R_frames2ref(:,:,i);
    pos = R_pos_ref(:,i);
    x = R_sz(i)*R*[1;0;0];
    y = R_sz(i)*R*[0;1;0];
    z = R_sz(i)*R*[0;0;1];
    quiver3(pos(1),pos(2),pos(3),x(1),x(2),x(3),'Color',c,'LineWidth',R_lw(i))
    quiver3(pos(1),pos(2),pos(3),y(1),y(2),y(3),'Color',c,'LineWidth',R_lw(i))
    quiver3(pos(1),pos(2),pos(3),z(1),z(2),z(3),'Color',c,'LineWidth',R_lw(i))
    x = x + pos;
    y = y + pos;
    z = z + pos;
    text(x(1),x(2),x(3),['$X_{',R_lg{i},'}$'],'Color',c,'Interpreter','latex')
    text(y(1),y(2),y(3),['$Y_{',R_lg{i},'}$'],'Color',c,'Interpreter','latex')
    text(z(1),z(2),z(3),['$Z_{',R_lg{i},'}$'],'Color',c,'Interpreter','latex')
end

for j = 1:m
    c = col(j+i,:);
    v = v_ref(:,j);
    pos = v_pos_ref(:,j);
    quiver3(pos(1),pos(2),pos(3),v(1),v(2),v(3),'Color',c,'LineWidth',v_lw(j))
    v = v + pos;
    text(v(1),v(2),v(3),['$',v_lg{j},'$'],'Color',c,'Interpreter','latex')
end


end