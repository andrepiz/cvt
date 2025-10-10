function ecFrame = accumulateFrames(ecQuant, ixFrames)
% ACCUMULATEFRAMES Accumulates slices of a 3D matrix using start or [start, end] indices.
%
%   ecFrame = accumulateFrames(ecQuant, ixFrames)
%
%   Inputs:
%       ecQuant  - MxNxT matrix
%       ixFrames - 
%         - Sx1 vector of start indices (old behavior)
%         - OR Sx2 matrix of [startIdx, endIdx] per frame (new behavior)
%
%   Output:
%       ecFrame  - MxNx(S or S-1) matrix with each slice being sum of frames
%                  over specified intervals

    [M, N, T] = size(ecQuant);

    % --- Single-column input: treat as start indices ---
    if isvector(ixFrames) || size(ixFrames, 2) == 1
        ixStart = ixFrames(:);
        S = length(ixStart);

        if S < 2
            error('ixFrames must contain at least two start indices.');
        end

        ecFrame = zeros(M, N, S);

        for i = 1:S
            startIdx = ixStart(i);
            if i == S
                endIdx = T;
            else
                endIdx = ixStart(i + 1) - 1;
            end

            if startIdx > T
                warning('Start index %d exceeds number of frames (%d). Skipping.', startIdx, T);
                continue;
            end
            startIdx = max(startIdx, 1);
            endIdx   = min(endIdx, T);

            if endIdx < startIdx
                warning('Start index %d is after end index %d. Skipping.', startIdx, endIdx);
                continue;
            end

            ecFrame(:, :, i) = sum(ecQuant(:, :, startIdx:endIdx), 3);
        end

    % --- Two-column input: explicit [start, end] pairs ---
    elseif size(ixFrames, 2) == 2
        S = size(ixFrames, 1);
        ecFrame = zeros(M, N, S);

        for i = 1:S
            startIdx = ixFrames(i, 1);
            endIdx   = ixFrames(i, 2);

            if startIdx > T
                warning('Start index %d exceeds number of frames (%d). Skipping.', startIdx, T);
                continue;
            end
            startIdx = max(startIdx, 1);
            endIdx   = min(endIdx, T);

            if endIdx < startIdx
                warning('End index %d is before start index %d. Skipping.', endIdx, startIdx);
                continue;
            end

            ecFrame(:, :, i) = sum(ecQuant(:, :, startIdx:endIdx), 3);
        end
    else
        error('ixFrames must be either a column vector or an Sx2 matrix of [start, end] indices.');
    end
end
