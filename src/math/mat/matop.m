function mat_out = matop(mat_in, op, varargin)
% Apply transformation to a matrix of 4 types:
% 'merge' ['sum'/'mean']: apply specified operation along the third dimension
% 'scale' [k]: multiply every value by k
% 'domain' [l, u]: re-map to a new domain with l and u as bounds
% 'gamma' [g]: multiply to the power of g
% 'mean' [m]: scale such that mean becomes m
% 'shift' [s]: shift of a quantity s
% If varargin is empty, returns the original matrix

if iscell(op) && length(varargin)==length(op)
    % apply transformations one consequent to other
    mat_out = mat_in;
    for ix = 1:length(op)
        mat_out = matop(mat_out, op{ix}, varargin{ix});
    end
    return
end

if isempty(varargin{1})
    mat_out = mat_in;
    return
end

switch op
    case 'merge'
        switch varargin{1}
            case 'sum'
                mat_out = sum(mat_in, 3);
            case 'mean'
                mat_out = mean(mat_in, 3);
            otherwise
                error('Merge transformation not recognized')
        end
    case 'scale'
        if isnumeric(varargin{1})
            mat_out = varargin{1}.*mat_in;
        else
            error('Scale factor is not a number')
        end
    case 'shift'
        if isnumeric(varargin{1})
            mat_out = varargin{1}+mat_in;
        else
            error('Shift factor is not a number')
        end
    case 'domain'
        m_in = min(mat_in, [],'all');
        M_in = max(mat_in, [],'all');
        if ~isnumeric(varargin{1})
            error('Domain parameter must be a vector of 2 elements')
        end
        m_out = varargin{1}(1);
        M_out = varargin{1}(2);
        mat_out = (mat_in - m_in)/(M_in - m_in)*(M_out - m_out) + m_out;
    case 'gamma'
        mat_out = mat_in.^(varargin{1});
    case 'mean'
        mat_out = mat_in*varargin{1}/mean(mat_in,'all');
    otherwise
        error('Operation not recognized')
end
