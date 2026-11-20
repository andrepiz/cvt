function fh = plotScene(ut, varargin)
%plotScene Visualize the scene
%
% Required:
%   ut                          1xN or Nx1 time vector
%   pos1, pos2, ..., posM       each 3xN position array in Moon-centered frame
%
% Optional name-value pairs:
%   'Pause'        scalar pause duration in seconds, default 0.1
%   'Colors'       Mx3 numeric RGB array, or cell array of color specs
%   'Names'        1xM cell array of names for legend
%   'Scale'        quiver length scale, default 20000e3
%   'MarkerSizes'  1xM vector of marker sizes
%
% Example:
%   fh = visualizeScene(ut, ...
%       pos_moon2sun_IAU, pos_moon2cam_IAU, pos_moon2earth_IAU, ...
%       'Pause', 0.05, ...
%       'Colors', {[1 0 0], [0 0.6 0], [0 0 1]}, ...
%       'Names', {'Sun','S/C','Earth'});

    % Find where name-value pairs start
    firstNV = find(cellfun(@(x) ischar(x) || isstring(x), varargin), 1, 'first');
    if isempty(firstNV)
        posArgs = varargin;
        nvArgs = {};
    else
        posArgs = varargin(1:firstNV-1);
        nvArgs = varargin(firstNV:end);
    end

    nBodies = numel(posArgs);
    if nBodies == 0
        error('At least one body position array must be provided.');
    end

    % Validate position inputs
    N = numel(ut);
    pos = cell(1, nBodies);
    for i = 1:nBodies
        pos{i} = posArgs{i};
        if ~isnumeric(pos{i}) || size(pos{i},1) ~= 3 || size(pos{i},2) ~= N
            error('Each position input must be a 3xN numeric array matching numel(ut).');
        end
    end

    % Defaults
    pauseTime = 0.1;
    scale = 20000e3;
    names = arrayfun(@(k) sprintf('Body %d', k), 1:nBodies, 'UniformOutput', false);
    markerSizes = 60 * ones(1, nBodies);

    defaultColors = lines(max(nBodies, 7));
    colors = defaultColors(1:nBodies, :);

    % Parse name-value pairs
    p = inputParser;
    addParameter(p, 'Pause', pauseTime, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Colors', colors);
    addParameter(p, 'Names', names, @(x) iscell(x) && numel(x) == nBodies);
    addParameter(p, 'Scale', scale, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'MarkerSizes', markerSizes, @(x) isnumeric(x) && numel(x) == nBodies);
    parse(p, nvArgs{:});

    pauseTime = p.Results.Pause;
    scale = p.Results.Scale;
    names = p.Results.Names;
    markerSizes = p.Results.MarkerSizes;
    colorsIn = p.Results.Colors;

    % Normalize color input
    if iscell(colorsIn)
        if numel(colorsIn) ~= nBodies
            error('If Colors is a cell array, it must contain one entry per body.');
        end
        colorsCell = colorsIn;
    elseif isnumeric(colorsIn) && isequal(size(colorsIn), [nBodies, 3])
        colorsCell = mat2cell(colorsIn, ones(1,nBodies), 3);
    else
        error('Colors must be either an Mx3 RGB array or a 1xM cell array of MATLAB color specs.');
    end

    % Precompute directions
    dirs = cell(1, nBodies);
    for i = 1:nBodies
        dirs{i} = vecnormalize(pos{i});
    end

    % Figure
    fh = figure( ...
        'Name', 'Scene', ...
        'Units', 'normalized', ...
        'Position', [0.1 0.1 0.6 0.5]);

    grid on
    hold on
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')

    scatter3(0, 0, 0, 80, 'k', 'filled', 'DisplayName', 'Origin');

    % Plot trajectories
    for i = 1:nBodies
        plot3(pos{i}(1,:), pos{i}(2,:), pos{i}(3,:), ...
            '-', 'Color', colorsCell{i}, ...
            'DisplayName', names{i});
    end
    legend show

    % Animated handles
    s = gobjects(1, nBodies);
    q = gobjects(1, nBodies);

    for ix = 1:N
        if ix > 1
            delete(s(ishandle(s)))
            delete(q(ishandle(q)))
        end

        for i = 1:nBodies
            s(i) = scatter3( ...
                pos{i}(1,ix), pos{i}(2,ix), pos{i}(3,ix), ...
                markerSizes(i), colorsCell{i}, 'filled', ...
                'HandleVisibility', 'off');

            q(i) = quiver3( ...
                0, 0, 0, ...
                scale * dirs{i}(1,ix), ...
                scale * dirs{i}(2,ix), ...
                scale * dirs{i}(3,ix), ...
                'Color', colorsCell{i}, ...
                'LineWidth', 1.5, ...
                'MaxHeadSize', 0.5, ...
                'HandleVisibility', 'off');
        end

        title(string(ut(ix)))
        drawnow

        if pauseTime > 0
            pause(pauseTime)
        end
    end
end