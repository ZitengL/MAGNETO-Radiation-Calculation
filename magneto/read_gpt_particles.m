function tr = read_gpt_particles(gdf_filepath, ids, skip_fE)
%% read_gpt_particles
% Reads the trajectories of selected particles from a (large) GPT tout output
% file without loading the whole file: the file is streamed block by block and
% only the requested particles (by GPT particle ID) and the needed arrays are kept.
% Use it to process a multi-particle run in batches of particles.
%
% Syntax:
%   tr = read_gpt_particles(gdf_filepath, ids)
%
% Inputs:
%   gdf_filepath - GPT output written with tout (one time group per output time)
%   ids          - GPT particle IDs to keep (vector)
%   skip_fE      - (optional) true: do not keep fEx, fEy, fEz (returned as zeros of
%                  size 1 x Nt); saves memory when there is no electric field. Default false.
%
% Output (struct): Nt = number of time groups, Np = numel(ids)
%   tr.t                - time [s], 1 x Nt
%   tr.ID               - particle IDs, Np x 1 (same order as ids)
%   tr.x, tr.y, tr.z    - positions [m], Np x Nt
%   tr.Bx, tr.By, tr.Bz - velocity / c, Np x Nt
%   tr.G                - Lorentz factor, Np x Nt
%   tr.fEx ... tr.fBz   - fields at the particles [V/m], [T], Np x Nt
%   tr.q, tr.m, tr.nmacro - charge [C], mass [kg], macro size, Np x 1 (from the first group)
% A particle missing from a group (lost) is NaN there.

ids = ids(:);
if nargin < 3; skip_fE = false; end
want = {'x', 'y', 'z', 'Bx', 'By', 'Bz', 'G', 'fEx', 'fEy', 'fEz', 'fBx', 'fBy', 'fBz'};
if skip_fE; want = setdiff(want, {'fEx', 'fEy', 'fEz'}, 'stable'); end
once = {'q', 'm', 'nmacro'};

GDFNAMELEN = 16; GDFID = 94325877;
t_dir = 256; t_edir = 512; t_param = 1024; t_data = 2048;
t_ascii = 1; t_s32 = 2; t_dbl = 3; t_nul = 16;

fid = fopen(gdf_filepath, 'r', 'l');
if fid < 0; error('Cannot open %s', gdf_filepath); end
cleaner = onCleanup(@() fclose(fid));
if fread(fid, 1, '*uint32') ~= GDFID; error('%s is not a GDF file', gdf_filepath); end
fseek(fid, 4 + 2*GDFNAMELEN + 8, 'cof');                % rest of the main header

% grow storage in chunks of time groups
chunk = 2000; cap = chunk; ng = 0;
for k = 1:numel(want); tr.(want{k}) = nan(numel(ids), cap); end
t = nan(1, cap);
for k = 1:numel(once); tr.(once{k}) = nan(numel(ids), 1); end
cur_time = NaN; idx = []; got_once = false; in_group = false;

while true
    name = fread(fid, GDFNAMELEN, 'uint8=>char')';          % not '*char': UTF-8 decoding runs out of memory on files of > ~20 GB
    if numel(name) < GDFNAMELEN; break; end              % end of file
    name = name(1:find([name 0] == 0, 1) - 1);
    block_type = fread(fid, 1, 'uint32');
    block_size = fread(fid, 1, 'uint32');
    dtype = bitand(block_type, 255);

    if bitand(block_type, t_dir); in_group = true; end
    if bitand(block_type, t_param)
        switch dtype
            case t_dbl
                v = fread(fid, 1, 'double');
                if strcmp(name, 'time'); cur_time = v; end
            case t_s32;   fread(fid, 1, 'int32');
            case t_ascii; fread(fid, block_size, 'uint8=>char');
            case t_nul
            otherwise; error('Unknown parameter type at byte %d', ftell(fid));
        end
    end
    if bitand(block_type, t_data)
        if dtype ~= t_dbl; error('Unknown array type at byte %d', ftell(fid)); end
        if strcmp(name, 'ID')
            v = fread(fid, block_size/8, 'double');
            [~, idx] = ismember(ids, v);                 % positions of the wanted particles in this group
        elseif any(strcmp(name, want)) || (~got_once && any(strcmp(name, once)))
            v = fread(fid, block_size/8, 'double');
            arrays_pending.(name) = v; %#ok<STRNU>
        else
            fseek(fid, block_size, 'cof');
        end
    end
    if bitand(block_type, t_edir) && in_group && exist('arrays_pending', 'var')
        % end of a time group: store the wanted particles
        if isempty(idx); error('No ID array in the group at t = %g', cur_time); end
        ng = ng + 1;
        if ng > cap
            cap = cap + chunk; t(cap) = NaN;
            for k = 1:numel(want); tr.(want{k})(:, cap) = NaN; end
        end
        t(ng) = cur_time;
        ok = idx > 0;
        for k = 1:numel(want)
            if isfield(arrays_pending, want{k}); tr.(want{k})(ok, ng) = arrays_pending.(want{k})(idx(ok)); end
        end
        if ~got_once
            for k = 1:numel(once)
                if isfield(arrays_pending, once{k}); tr.(once{k})(ok) = arrays_pending.(once{k})(idx(ok)); end
            end
            got_once = true;
        end
        clear arrays_pending
        idx = []; in_group = false;
    end
end
tr.t = t(1:ng);
for k = 1:numel(want); tr.(want{k}) = tr.(want{k})(:, 1:ng); end
if skip_fE; tr.fEx = zeros(1, ng); tr.fEy = zeros(1, ng); tr.fEz = zeros(1, ng); end
tr.ID = ids;
tr.source = gdf_filepath;
end
