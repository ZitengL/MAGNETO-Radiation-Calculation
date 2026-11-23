%% compute_bunch_vs_z
% Bunch length and longitudinal form factor along the beamline, from the GPT particle output
% (read only, streamed one time group at a time; only the z array is read):
%   at each GPT output time t:  <z>,  sigma_z = std(z_j),  |F|^2 = |(1/N) sum_j exp(i k z_j)|^2,  k = 2 pi / 70 nm
% |F|^2 is the form factor of the actual 2000-electron sample (not a Gaussian formula); its floor for an
% incoherent sample is ~1/N = 5e-4. Positions are taken at a fixed time t (GPT tout output).
% Cases: sd3_en0 (sigma_delta = 3 %) and sd0_en12 (epsilon_n = 12 um); sd0_en0 (cold) as reference.
% Output: data/spread_emittance/bunch_vs_z_<case>.txt and .mat
% Run: matlab -batch "compute_bunch_vs_z"

this_dir = fileparts(mfilename('fullpath'));
gpt_dir  = fullfile(this_dir, '..', '..', 'gpt_sims', 'multi_electron', 'E240MeV_und15p_N2000_spread_emittance');
out_dir  = fullfile(this_dir, 'data', 'spread_emittance');
lam_F = 70e-9; k = 2*pi/lam_F;
if ~exist('cases', 'var'); cases = {'sd3_en0', 'sd0_en12', 'sd0_en0'}; end   % override: matlab -batch "cases={...}; compute_bunch_vs_z"

addpath(fullfile(this_dir, '..', '..', 'magneto'));                       % load_gdf.m (for the small stats.gdf check)
for c = 1:numel(cases)
    f = fullfile(gpt_dir, cases{c}, 'result.gdf');
    tic; [t, zm, sz, F2, np] = stream_z(f, k);
    fprintf('%s: %d time groups, N = %d, %.0f s\n', cases{c}, numel(t), np(1), toc);
    % check against GPT's own statistics (stats.gdf, same output times)
    g = load_gdf(fullfile(gpt_dir, cases{c}, 'stats.gdf')); s = g(end).d;
    sz_gpt = interp1(s.time(:), s.stdz(:), t, 'linear', NaN);
    ok = ~isnan(sz_gpt);
    fprintf('   max |sigma_z - GPT stdz| / sigma_z = %.2e (on %d times)\n', max(abs(sz(ok) - sz_gpt(ok)) ./ sz(ok)), nnz(ok));
    M = [t(:), zm(:), sz(:), F2(:), np(:)];
    save(fullfile(out_dir, ['bunch_vs_z_' cases{c} '.mat']), 'M', 'lam_F');
    fid = fopen(fullfile(out_dir, ['bunch_vs_z_' cases{c} '.txt']), 'w');
    fprintf(fid, '%% %s, GPT output %s (read only). At each output time: <z>, sigma_z = std(z_j), |F|^2 = |(1/N) sum exp(i k z_j)|^2, k = 2 pi/%.0f nm\n', ...
        cases{c}, f, lam_F*1e9);
    fprintf(fid, '%% Columns: 1 t_s  2 mean_z_m  3 sigma_z_m  4 F2_70nm  5 N\n');
    fprintf(fid, '%.17e %.17e %.17e %.17e %d\n', M');
    fclose(fid);
end

function [t, zm, sz, F2, np] = stream_z(f, k)
GDFNAMELEN = 16; GDFID = 94325877;
t_dir = 256; t_edir = 512; t_param = 1024; t_data = 2048;
t_ascii = 1; t_s32 = 2; t_dbl = 3; t_nul = 16;
fid = fopen(f, 'r', 'l'); if fid < 0; error('Cannot open %s', f); end
cleaner = onCleanup(@() fclose(fid));
if fread(fid, 1, '*uint32') ~= GDFID; error('%s is not a GDF file', f); end
fseek(fid, 4 + 2*GDFNAMELEN + 8, 'cof');
cap = 50000; t = nan(1, cap); zm = t; sz = t; F2 = t; np = t; ng = 0;
cur_time = NaN; z = []; in_group = false;
while true
    name = fread(fid, GDFNAMELEN, 'uint8=>char')';
    if numel(name) < GDFNAMELEN; break; end
    name = name(1:find([name 0] == 0, 1) - 1);
    block_type = fread(fid, 1, 'uint32'); block_size = fread(fid, 1, 'uint32');
    dtype = bitand(block_type, 255);
    if bitand(block_type, t_dir); in_group = true; end
    if bitand(block_type, t_param)
        switch dtype
            case t_dbl; v = fread(fid, 1, 'double'); if strcmp(name, 'time'); cur_time = v; end
            case t_s32; fread(fid, 1, 'int32');
            case t_ascii; fread(fid, block_size, 'uint8=>char');
            case t_nul
            otherwise; error('Unknown parameter type at byte %d', ftell(fid));
        end
    end
    if bitand(block_type, t_data)
        if dtype ~= t_dbl; error('Unknown array type at byte %d', ftell(fid)); end
        if strcmp(name, 'z'); z = fread(fid, block_size/8, 'double'); else; fseek(fid, block_size, 'cof'); end
    end
    if bitand(block_type, t_edir) && in_group && ~isempty(z)
        ng = ng + 1;
        if ng > cap; cap = cap + 50000; t(cap) = NaN; zm(cap) = NaN; sz(cap) = NaN; F2(cap) = NaN; np(cap) = NaN; end
        mz = mean(z); dz = z - mz;                           % subtract the mean before the phase: full precision
        t(ng) = cur_time; zm(ng) = mz; sz(ng) = std(z, 1); np(ng) = numel(z);
        F2(ng) = abs(mean(exp(1i * k * dz)))^2;
        z = []; in_group = false;
    end
end
t = t(1:ng); zm = zm(1:ng); sz = sz(1:ng); F2 = F2(1:ng); np = np(1:ng);
end
