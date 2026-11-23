function [gdf, info] = load_gdf(filename, arrays_to_load, elements_to_load)
%LOAD_GDF Load a binary GPT GDF file into a MATLAB structure array.
%
% This function was copied from the local load_gdf function in
% Brookhaven undulator/ATF_IFEL.m so that MAGNETO can call it.

if nargin == 1
    arrays_to_load   = []; % load all arrays
    elements_to_load = []; % load all elements
    last_element_to_load = -1; % load all elements
elseif nargin == 2
    elements_to_load = []; % load all elements
    last_element_to_load = -1; % load all elements
else
    last_element_to_load = max(elements_to_load); % stop when this element is loaded
end

% Constants
GDFNAMELEN = 16;                 % Length of the ascii-names
GDFID  = 94325877;

% Data types
t_ascii  = hex2dec('0001');      % Char array
t_s32    = hex2dec('0002');      % Signed long
t_dbl    = hex2dec('0003');      % Double
t_nul    = hex2dec('0010');      % No data

% Block types
t_dir    = hex2dec('0100');      % Directory entry start
t_edir   = hex2dec('0200');      % Directory entry end
t_param  = hex2dec('0400');      % Parameter
t_data   = hex2dec('0800');      % Data array

% Open file
FileID = fopen(filename, 'r', 'l', 'windows-1252');
if FileID < 0
    error('Error in ''load_gdf'': did not find file ''%s''', filename);
end

% Read gdfmainheader
ID = fread(FileID, 1, '*uint32');
if ID ~= GDFID
    fclose(FileID);
    error('Error in ''load_gdf'': this is not a gdf file');
end
info = struct;
info.ID = ID;
info.cretime  = fread(FileID, 1, '*uint32');
info.creator  = fread(FileID, GDFNAMELEN, '*char')';
info.creator  = info.creator(1:find(info.creator == 0, 1)-1);
info.destin   = fread(FileID, GDFNAMELEN, '*char')';
info.destin   = info.destin(1:find(info.destin == 0, 1)-1);
info.gdfmaj   = fread(FileID, 1, '*uint8');
info.gdfmin   = fread(FileID, 1, '*uint8');
info.cremaj   = fread(FileID, 1, '*uint8');
info.cremin   = fread(FileID, 1, '*uint8');
info.desmaj   = fread(FileID, 1, '*uint8');
info.desmin   = fread(FileID, 1, '*uint8');
dummy         = fread(FileID, 2, '*uint8'); %#ok<NASGU>

% Read gdf data blocks
gdf    = struct;
params = struct;
params.level1 = [];
arrays = struct;
current_element = 1;
elements_stored = 0;
level = 1;

while ~feof(FileID)
    % Read gdf block header
    name = fread(FileID, GDFNAMELEN, '*char')';
    name = name(1:find(name == 0, 1)-1);

    % Get rid of not allowed characters etc.
    if isempty(name)
        name = 'empty';
    elseif name(1) == '@'
        name = 'x0x40logo';
    end

    block_type = fread(FileID, 1, '*uint32');
    data_type  = bitand(block_type, 255);
    block_size = fread(FileID, 1, '*uint32');

    start_dir     = bitand(block_type, t_dir) > 0;
    end_dir       = bitand(block_type, t_edir) > 0;
    data_is_param = bitand(block_type, t_param) > 0;
    data_is_array = bitand(block_type, t_data) > 0;

    if start_dir
        level = level + 1;
        params.(sprintf('level%d', level)) = ...
            params.(sprintf('level%d', level-1));
    end

    if end_dir
        if isempty(fieldnames(arrays))
            if isempty(elements_to_load)
                if ~iscell(arrays_to_load) && strcmp(arrays_to_load, ' ')
                    elements_stored = elements_stored + 1;
                elseif elements_stored == 0
                    elements_stored = 1;
                end
                gdf(elements_stored).p = params.(sprintf('level%d', level));
            end
        else
            if isempty(elements_to_load) || any(elements_to_load == current_element)
                elements_stored = elements_stored + 1;
                gdf(elements_stored).d = arrays;
                gdf(elements_stored).p = params.(sprintf('level%d', level));
                if current_element == last_element_to_load
                    break
                end
            end
            arrays = struct;
            current_element = current_element + 1;
        end
        level = level - 1;
    end

    % Read a scanned parameter, time output, or screen position.
    if data_is_param
        switch data_type
            case t_dbl
                value = fread(FileID, 1, 'double');
                params.(sprintf('level%d', level)).(name) = value;
            case t_nul
                % No data present.
            case t_ascii
                value = fread(FileID, double(block_size), '*char')';
                params.(sprintf('level%d', level)).(name) = value;
            case t_s32
                value = fread(FileID, 1, 'int32')';
                params.(sprintf('level%d', level)).(name) = value;
            otherwise
                fclose(FileID);
                error(['Error in ''load_gdf'': unknown datatype of value\n' ...
                    'Byte position in file: %d\n'], ftell(FileID));
        end
    end

    % Read particle data arrays.
    if data_is_array
        switch data_type
            case t_dbl
                if mod(block_size, 8) ~= 0
                    fclose(FileID);
                    error('Error in ''load_gdf'': wrong size for double array');
                end
                N = double(block_size/8);
                value = fread(FileID, N, 'double');

                if isempty(arrays_to_load) || any(strcmp(name, arrays_to_load))
                    try
                        arrays.(name) = value;
                    catch
                        name = matlab.lang.makeValidName(name);
                        arrays.(name) = value;
                    end
                end
            otherwise
                fclose(FileID);
                error(['Error in ''load_gdf'': unknown datatype of array\n' ...
                    'Byte position in file: %d\n'], ftell(FileID));
        end
    end
end

% params contains numderivs and cpu for the last element.
if isempty(elements_to_load) && isfield(params, sprintf('level%d', level)) && ...
        isfield(params.(sprintf('level%d', level)), 'cputime')
    if elements_stored == 0
        elements_stored = 1;
    end
    level_str = sprintf('level%d', level);
    field_names = fieldnames(params.(level_str));
    for i = 1:length(field_names)
        gdf(elements_stored).p.(field_names{i}) = ...
            params.(level_str).(field_names{i});
    end
end

% Some gdf files do not use start/end directory identifiers.
if ~isempty(fieldnames(arrays)) && isempty(elements_to_load)
    elements_stored = elements_stored + 1;
    gdf(elements_stored).d = arrays;
    gdf(elements_stored).p = params.(sprintf('level%d', level));
end

fclose(FileID);

end
