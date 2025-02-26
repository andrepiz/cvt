function p = plotBooleanCube(position, size, bool)
    % PLOTBOOLEANCUBE Plots a semi-transparent red/green cube.
    % 
    % position - 3-element vector [x, y, z] specifying the cube's bottom-left-front corner.
    % size - Scalar or 3-element vector [sx, sy, sz] specifying the cube dimensions.
    
    if numel(size) == 1
        size = [size, size, size]; % Make it a cube if scalar is provided
    end
    
    % Define cube vertices based on position and size
    x = [0 1 1 0 0 1 1 0] * size(1) + position(1);
    y = [0 0 1 1 0 0 1 1] * size(2) + position(2);
    z = [0 0 0 0 1 1 1 1] * size(3) + position(3);

    % Define faces
    faces = [1 2 3 4;  % Bottom
             5 6 7 8;  % Top
             1 2 6 5;  % Front
             2 3 7 6;  % Right
             3 4 8 7;  % Back
             4 1 5 8]; % Left

    % Plot the cube
    if bool
        col = 'g';
    else
        col = 'r';
    end
    p = patch('Vertices', [x(:), y(:), z(:)], 'Faces', faces, ...
          'FaceColor', col, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'LineWidth', 1.5);
end