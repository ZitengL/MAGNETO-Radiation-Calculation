function run_direct_sigz(sz_nm, chunks, gpt_dir_name, nthreads)
%% run_direct_sigz
% User (2026-11-06): full 5e4-particle sigma_z scan computed DIRECTLY (a new GPT run for every sigma_z, no time shift,
% no expected value), at x = 0, y = -2.0 mm, z = 5.9 m, lambda 50-1000 nm; raw |sum E|^2 shown with the macro-particle
% noise floor. Stop where the direct result can no longer resolve sigma_z (user: no need to go to 1000 nm).
% Grid chosen from the expected-value result vs the macro-particle noise floor (sum |E_j|^2 = nmacro x incoherent):
%   D1 (50-280 nm): coherent = noise floor near 70 nm; D2 (300-1000 nm): near 150 nm  ->  20:10:150 nm and 200 nm.
%   sigma_z = 10 nm is the GPT run already made in run_large_sigz.m (Er at 10 nm has no phase factor = direct).
% Everything else (beam, random numbers, GPT input, dz = 120 um, undulator/dipole split, lambda grid) as run_large_sigz.m.
% Stored per sigma_z: results/direct_sigz/direct_sz<nnnn>nm.mat, struct D:
%   Sd(chunk, part, comp, w) = sum_j E_j   (part 1 undulator z < 0.63 m, 2 dipole; comp x y z; fields carry q*nmacro)
%   I (chunk, part, comp, w) = sum_j |E_j|^2 (part 1 und, 2 dip, 3 und + dip)
% Resumable (per sigma_z and chunk). Progress: results/direct_sigz/progress.log
% Usage: run_direct_sigz([20:10:150 200], 1:25)
% Two MATLAB at once (2026-11-07, user asleep): each needs its OWN GPT folder (chunk_beam/chunk_result are fixed names)
% and a DISJOINT sigma_z list (one .mat per sigma_z, no locking), e.g.
%   run_direct_sigz([30 50 70], 1:25)                                   (gpt_sims/multi_electron/E7GeV_und_dipole)
%   run_direct_sigz([100 150 200], 1:25, 'E7GeV_und_dipole_w2', 8)     (identical copy of multi_electron.in + run.bat)

if nargin < 3 || isempty(gpt_dir_name); gpt_dir_name = 'E7GeV_und_dipole'; end
if nargin < 4 || isempty(nthreads); nthreads = 8; end
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
% (read_gpt_particles and particle_traj are in magneto/)
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', gpt_dir_name);
out_dir = fullfile(this_dir, 'results', 'direct_sigz'); if ~exist(out_dir, 'dir'); mkdir(out_dir); end
maxNumCompThreads(nthreads);

c0 = 299792458; dz_um = 120; zsplit = 0.63;
obs = [0 -2e-3 5.9];
lam = logspace(log10(50), log10(1000), 600)*1e-9; om = 2*pi*c0./lam; nw = numel(lam);
Ntot = 5e4; Nc = 2000; nch = Ntot/Nc;
S = load(fullfile(this_dir, '..', 'E7GeV_gap6und_sigz_sweep', 'results', 'standard_sample.mat'));
S.X(:, 1:2) = -S.X(:, 1:2); S.U(:, 1:2) = -S.U(:, 1:2);           % rotation by 180 deg about z (as run_large_sigz)
logf = fullfile(out_dir, 'progress.log');

for sz = sz_nm(:)'
    out = fullfile(out_dir, sprintf('direct_sz%04dnm.mat', sz));
    if exist(out, 'file'); D = load(out); D = D.D;
    else
        D = struct('sigma_z', sz*1e-9, 'obs', obs, 'lambda', lam, 'omega', om, 'dz_um', dz_um, 'nmacro', S.nmacro, ...
            'Sd', zeros(nch, 2, 3, nw), 'I', zeros(nch, 3, 3, nw), 'done', false(1, nch), 'tchunk', nan(1, nch), ...
            'note', 'direct GPT run at sigma_z; Sd = sum_j E_j, I = sum_j |E_j|^2; part 1 und, 2 dip (3 = und + dip for I)');
    end
    for ic = chunks
        if D.done(ic); continue; end
        tc = tic; ids = (ic-1)*Nc + (1:Nc);
        D0 = write_chunk(gdir, S, ids, sz*1e-9, 0.4);
        gdf = fullfile(gdir, 'chunk_result.gdf'); if exist(gdf, 'file'); delete(gdf); end
        [st, msg] = system(sprintf('cmd /c ""%s" %d"', fullfile(gdir, 'run.bat'), dz_um));
        if st ~= 0 || ~exist(gdf, 'file'); error('GPT failed (sz %d nm, chunk %d): %s', sz, ic, msg); end
        tg = toc(tc);
        Sd = zeros(2, 3, nw); I = zeros(3, 3, nw);
        for b = 1:2
            loc = (b-1)*1000 + (1:1000);
            tr = read_gpt_particles(gdf, loc, true);
            zp = tr.z(:,1); xp = D0(loc,1) + D0(loc,4)./D0(loc,6).*(zp - D0(loc,3)); yp = D0(loc,2) + D0(loc,5)./D0(loc,6).*(zp - D0(loc,3));
            if max(abs([tr.x(:,1) - xp; tr.y(:,1) - yp])) > 1e-9; error('ID / initial distribution mismatch'); end
            for k = 1:numel(loc)
                traj = particle_traj(tr, k);
                Eu = emission_time_spectrum(sub(traj, traj.r(3,:) < zsplit), obs, lam, 'fB', norm(obs)/c0).Ew;
                Ed = emission_time_spectrum(sub(traj, traj.r(3,:) >= zsplit), obs, lam, 'fB', norm(obs)/c0).Ew;
                Sd(1,:,:) = squeeze(Sd(1,:,:)) + Eu;  Sd(2,:,:) = squeeze(Sd(2,:,:)) + Ed;
                I(1,:,:) = squeeze(I(1,:,:)) + abs(Eu).^2;
                I(2,:,:) = squeeze(I(2,:,:)) + abs(Ed).^2;
                I(3,:,:) = squeeze(I(3,:,:)) + abs(Eu + Ed).^2;
            end
            clear tr
        end
        delete(gdf);
        D.Sd(ic,:,:,:) = Sd; D.I(ic,:,:,:) = I; D.done(ic) = true; D.tchunk(ic) = toc(tc);
        save(out, 'D');
        fid = fopen(logf, 'a');
        if fid > 0       % two workers may append at the same moment; never stop the run over a log line
            fprintf(fid, 'sz %4d nm, chunk %2d/%d: GPT %.0f s, total %.0f s  [%s]  (%s)\n', sz, ic, nch, tg, toc(tc), ...
                datestr(now, 'mm-dd HH:MM'), gpt_dir_name);
            fclose(fid);
        end
    end
end
end

function D0 = write_chunk(gdir, S, ids, sz, zw)
% identical to run_large_sigz.m / run_combined.m
X = S.X(ids,:); U = S.U(ids,:); n = numel(ids);
xw = S.sx*X(:,1); yw = S.sx*X(:,2); z = sz*X(:,3);
xp = S.sth*U(:,1); yp = S.sth*U(:,2); G = S.G0*(1 + S.sig_d*U(:,3));
GB = sqrt(G.^2 - 1); GBz = GB./sqrt(1 + xp.^2 + yp.^2);
D0 = [xw - xp*zw, yw - yp*zw, z, xp.*GBz, yp.*GBz, GBz, G, S.nmacro*ones(n, 1)];
fid = fopen(fullfile(gdir, 'chunk_beam.txt'), 'wt');
fprintf(fid, 'x y z GBx GBy GBz G nmacro\n'); fprintf(fid, '%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n', D0'); fclose(fid);
[st, msg] = system(sprintf('cd /d "%s" && "%%GPTBIN%%\\asci2gdf.exe" -o chunk_beam.gdf chunk_beam.txt', gdir));
if st ~= 0; error('asci2gdf failed: %s', msg); end
end

function t = sub(t, m)
% identical to run_large_sigz.m
k = find(m); kf = k(any(t.fB(:, k) ~= 0, 1)); k = max(k(1), kf(1)-1) : min(k(end), kf(end)+1);
for f = {'t', 'G'};             t.(f{1}) = t.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; t.(f{1}) = t.(f{1})(:, k); end
end
