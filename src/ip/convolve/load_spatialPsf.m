function [psfs, spacing, offset_horizontal, offset_vertical] = load_spatialPsf(psfs, spacing, offset_horizontal, offset_vertical, source)
%LOAD_SPATIALPSF Load a spatial PSF grid and resolve optional metadata.
%   PSFS is an Nv-by-Nu cell array of numeric intensity matrices or file
%   names (mixed cells are allowed). SOURCE selects the file reader:
%     'zemax' - Huygens PSF text export, including its header (default).
%     'csv'   - numeric comma-delimited data with decimal points, no header.
%   Nonempty metadata arguments always take precedence over file headers.
%   Empty SPACING is inferred in metres; empty offsets are inferred in
%   radians from Zemax's degree-valued field coordinates. Numeric-only/CSV
%   grids require explicit metadata. Mixed numeric/Zemax grids can infer
%   metadata from the available headers. One spacing applies to the grid.
%   A horizontal node belongs to each cell column, a vertical node to each
%   cell row. Inferred coordinates must form a rectangular monotonic grid.
%   Kernel arrays and coordinate signs are preserved: no rotation, flip or
%   normalization is applied here. Use the same right/down convention as
%   APPLY_SPATIALPSF and arrange filenames in their intended grid positions.

if nargin < 5 || isempty(source)
    source = 'zemax';
end
if nargin < 4, offset_vertical = []; end
if nargin < 3, offset_horizontal = []; end
if nargin < 2, spacing = []; end
assert(iscell(psfs) && ismatrix(psfs) && ~isempty(psfs), ...
    'load_spatialPsf:InvalidGrid', 'kernel must be a nonempty 2-D cell array.');
assert((ischar(source) && isrow(source)) || (isstring(source) && isscalar(source)), ...
    'load_spatialPsf:InvalidSource', 'source must be a character vector or string scalar.');
source = lower(strtrim(char(source)));
[nv, nu] = size(psfs);
header_spacing = nan(nv, nu);
header_x = nan(nv, nu);
header_y = nan(nv, nu);

for v = 1:nv
    for u = 1:nu
        value = psfs{v,u};
        if (ischar(value) && isrow(value)) || (isstring(value) && isscalar(value))
            filename = char(value);
            switch source
                case 'zemax'
                    [value, ds, x, y] = read_zemax_psf(filename, ...
                        isempty(spacing), isempty(offset_horizontal), isempty(offset_vertical));
                    header_spacing(v,u) = ds;
                    header_x(v,u) = x;
                    header_y(v,u) = y;
                case 'csv'
                    value = readmatrix(filename, 'FileType','text', ...
                        'Delimiter',',', 'DecimalSeparator','.');
                otherwise
                    error('load_spatialPsf:UnsupportedSource', ...
                        'Unsupported PSF source "%s". Use "zemax" or "csv".', source);
            end
        end
        validateattributes(value, {'numeric'}, ...
            {'real','finite','2d','nonempty','nonnegative'}, mfilename, 'kernel');
        assert(any(value(:)>0), 'load_spatialPsf:ZeroKernel', ...
            'PSF at row %d, column %d has zero total intensity.', v, u);
        psfs{v,u} = double(value);
    end
end

if isempty(spacing)
    available = header_spacing(isfinite(header_spacing));
    assert(~isempty(available), 'load_spatialPsf:MissingSpacing', ...
        'Provide spacing in metres, or supply a Zemax export containing data spacing.');
    assert_consistent(available, 'PSF sample spacing');
    spacing = available(1);
end
validateattributes(spacing, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'spacing');
if isempty(offset_horizontal)
    offset_horizontal = infer_nodes(header_x, 'horizontal');
end
if isempty(offset_vertical)
    offset_vertical = infer_nodes(header_y.', 'vertical');
end
validate_nodes(offset_horizontal, nu, 'offset_horizontal');
validate_nodes(offset_vertical, nv, 'offset_vertical');
end

function nodes = infer_nodes(metadata, direction)
nodes = zeros(1,size(metadata,2));
for i = 1:numel(nodes)
    available = metadata(isfinite(metadata(:,i)),i);
    assert(~isempty(available), 'load_spatialPsf:MissingOffsets', ...
        'Cannot infer %s node %d; provide offsets in radians.', direction, i);
    assert_consistent(available, [direction ' field angle']);
    nodes(i) = available(1);
end
end

function assert_consistent(values, description)
tolerance = max(1e-12, 1e-6*max(abs(values)));
assert(max(abs(values-values(1))) <= tolerance, ...
    'load_spatialPsf:InconsistentMetadata', ...
    'Inconsistent %s in the PSF grid. Check file positions or provide metadata explicitly.', description);
end

function validate_nodes(nodes, count, name)
validateattributes(nodes, {'numeric'}, ...
    {'real','finite','vector','numel',count}, mfilename, name);
assert(all(abs(nodes)<pi/2) && ...
    (count==1 || all(diff(nodes)>0) || all(diff(nodes)<0)), ...
    'load_spatialPsf:InvalidOffsets', ...
    '%s must be strictly monotonic field angles in radians between -pi/2 and pi/2.', name);
end

function [kernel, spacing, x, y] = read_zemax_psf(filename, need_spacing, need_x, need_y)
% Explicit BOM handling supports Zemax UTF-16 and ordinary UTF-8 exports.
fid = fopen(filename,'rb');
assert(fid >= 0, 'load_spatialPsf:FileOpen', 'Cannot open PSF file: %s', filename);
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid,Inf,'*uint8').';
if numel(bytes)>=2 && isequal(bytes(1:2),uint8([255 254]))
    txt = native2unicode(bytes(3:end),'UTF-16LE');
elseif numel(bytes)>=2 && isequal(bytes(1:2),uint8([254 255]))
    txt = native2unicode(bytes(3:end),'UTF-16BE');
else
    if numel(bytes)>=3 && isequal(bytes(1:3),uint8([239 187 191]))
        bytes = bytes(4:end);
    end
    txt = native2unicode(bytes,'UTF-8');
end
header = txt;
txt = strrep(txt,',','.');
number = '([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)';
spacing = NaN; x = NaN; y = NaN;
if need_spacing
    token = regexpi(txt, ['Data\s+spacing\s+is\s+' number '\s*([^\s.]+)'], 'tokens','once');
    assert(~isempty(token), 'load_spatialPsf:MissingSpacing', ...
        'No data spacing found in Zemax file: %s', filename);
    % Zemax Huygens headers report micrometres; also accept explicit units.
    unit = lower(token{2});
    switch unit
        case {'um', [char(181) 'm'], [char(956) 'm'], 'microns', 'micrometers', 'micrometres'}
            scale = 1e-6;
        case 'mm'
            scale = 1e-3;
        case 'nm'
            scale = 1e-9;
        case 'm'
            scale = 1;
        otherwise
            error('load_spatialPsf:SpacingUnit', ...
                'Unsupported spacing unit "%s" in %s; provide spacing explicitly.', unit, filename);
    end
    spacing = str2double(token{1})*scale;
    assert(isfinite(spacing) && spacing>0, 'load_spatialPsf:InvalidSpacing', ...
        'Invalid PSF sample spacing in %s.', filename);
end
if need_x || need_y
    % Example header: "0,4500 um at -1,5000, -1,5000 (deg)."
    field = extractBetween(header, ' at ', ' (deg)');
    assert(numel(field)==1, 'load_spatialPsf:MissingOffsets', ...
        'No field angles found in Zemax file: %s', filename);
    % Remove the comma-space between angles, then convert decimal commas.
    angles = sscanf(strrep(strrep(field{1}, ', ', ' '), ',', '.'), '%f');
    assert(numel(angles)==2, 'load_spatialPsf:MissingOffsets', ...
        'Expected two field angles in Zemax file: %s', filename);
    if need_x, x = deg2rad(angles(1)); end
    if need_y, y = deg2rad(angles(2)); end
end
shape = regexpi(txt,'Image\s+grid\s+size:\s*(\d+)\s+by\s+(\d+)','tokens','once');
assert(~isempty(shape), 'load_spatialPsf:MissingShape', ...
    'No image grid size found in Zemax file: %s', filename);
% Zemax reports horizontal samples first, then vertical samples.
nc = str2double(shape{1}); nr = str2double(shape{2});
assert(nr>0 && nc>0, 'load_spatialPsf:InvalidShape', 'Invalid grid size in %s.', filename);
data_start = regexpi(txt,'Values\s+are\s+relative\s+intensity\s*\.', 'end','once');
assert(~isempty(data_start), 'load_spatialPsf:MissingData', ...
    'Expected a relative-intensity Zemax PSF export: %s', filename);
values = sscanf(txt(data_start+1:end), '%f');
assert(numel(values)==nr*nc, 'load_spatialPsf:DataSize', ...
    'Expected %d matrix values in %s; found %d.', nr*nc, filename, numel(values));
kernel = reshape(values, nc, nr).';
end
