function noise = noise_smearing(ecr_in, direction, t_readout)
% Smearing is made of residual photons collected during the shift time from one row
% to the next of the signal at readout. The readout can be from up or down
% direction. This model assumes that the signal remains constant during
% readout, i.e., t_readout << t_exp
%
% Inputs:
%   ecr_in - Clean image electron count rate matrix (electrons/s), size [n_rows, n_cols]
%   direction - 'up' (bottom-to-top readout) or 'down' (top-to-bottom)
%   t_readout - Total readout time for full frame (s)
%
% Output:
%   noise - Smeared accumulation matrix (electrons), same size as ecr

[n_rows, ~] = size(ecr_in);

tshift   = t_readout / n_rows;
microexp = ecr_in * tshift;

switch direction
    case 'up'
        % smear from rows below → reverse cumulative sum
        noise = flipud(cumsum(flipud(microexp), 1)) - microexp;

    case 'down'
        % smear from rows above → forward cumulative sum
        noise = cumsum(microexp, 1) - microexp;

    otherwise
        error('direction must be ''up'' or ''down''');
end

end