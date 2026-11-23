function run_sweep(sz_nm_list, chunks, dz_um, tag)
%% run_sweep
% sigma_z sweep for the bunch-length idea (Neuman, Graves, O'Shea, PRSTAB 3, 030701 (2000)), SLAC gap-6 undulator,
% 7 GeV kinetic, 100 pC, point (b) = (0, 1.0 mm, 5.9 m) only. User's decisions (2026-11-04):
%   sigma_z = 10:5:50 nm; a GPT run for every sigma_z; 5e4 macro-particles; the SAME random numbers for every
%   sigma_z (only z is scaled); energy integrated over 50-300 nm only.
% Beam as gpt_sims/multi_electron/E7GeV_gap6und_N2000/make_distribution.m except N = 5e4 (nmacro = 1.248e4):
%   Gaussian, uncorrelated, cut at 4 sigma by redrawing; waist at z = 0.4 m with sigma_x = sigma_y = 20 um,
%   en = 30 um, sigma_delta = 2 % (dgamma/gamma); z at t = 0 with sigma_z; seeds 20261104 (x_w, y_w, z) and
%   20261105 (x', y', delta); straight-line drift from the start to the waist.
% Tracking in 25 chunks of 2000 particles (no space charge -> independent particles, identical to one GPT run);
% each chunk: write chunk_beam.gdf -> GPT (run.bat) -> MAGNETO -> delete chunk_result.gdf (disk: ~3 GB at a time).
% MAGNETO: emission_time_spectrum with t_ref = |obs|/c, fields carry q*nmacro; per chunk the sums
%   Ec(chunk, comp, omega) = sum_j E_j and Ic = sum_j |E_j|^2 are stored (25 independent subsets -> error bars).
% Wavelengths: 300 log points in 50-300 nm (band integral converged to 2e-5 against 2290 points, checked on final.mat).
% Runs at below-normal priority with at most 6 computational threads (the user works on the PC meanwhile).
% Resumable: results/<tag>sigz<NN>nm.mat is saved after every chunk.
%
% Syntax: run_sweep(sz_nm_list, chunks, dz_um, tag)   e.g. run_sweep(10:5:50, 1:25, 60, '')
%         tout check: run_sweep([10 50], 1, 30, 'dz30_')

this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
% (read_gpt_particles and particle_traj are in magneto/)                                      % read_gpt_particles, particle_traj
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_gap6und_sigz_sweep');
res_dir = fullfile(this_dir, 'results'); if ~exist(res_dir, 'dir'); mkdir(res_dir); end
maxNumCompThreads(6);
system(sprintf('powershell -NoProfile -Command "(Get-Process -Id %d).PriorityClass = ''BelowNormal''"', feature('getpid')));

c0 = 299792458; obs = [0 1e-3 5.9]; lam = logspace(log10(50), log10(300), 300)*1e-9;
Ntot = 5e4; Nc = 2000; nch = Ntot/Nc; zw = 0.4; zund = 0.4; zmap = 0.22381;
S = standard_sample(res_dir, Ntot);

for sz_nm = sz_nm_list
    out = fullfile(res_dir, sprintf('%ssigz%02dnm.mat', tag, sz_nm));
    if exist(out, 'file'); R = load(out); R = R.R;
    else
        R = struct('sigma_z', sz_nm*1e-9, 'obs', obs, 'lambda', lam, 'omega', 2*pi*c0./lam, 'dz_um', dz_um, ...
            'Ec', zeros(nch, 3, numel(lam)), 'Ic', zeros(nch, 3, numel(lam)), 'done', false(1, nch), 'xmax', 0, ...
            'nmacro', S.nmacro, 'tchunk', nan(1, nch), ...
            'note', 'Ec(chunk,comp,omega) = sum_j E_j (fields carry q*nmacro), Ic = sum_j |E_j|^2; MAGNETO = eps0 c |obs|^2/pi |sum_chunks Ec|^2');
    end
    for ic = chunks
        if R.done(ic); continue; end
        tc = tic; ids = (ic-1)*Nc + (1:Nc);
        D0 = write_chunk(gdir, S, ids, sz_nm*1e-9, zw);
        [st, msg] = system(sprintf('cmd /c ""%s" %d"', fullfile(gdir, 'run.bat'), dz_um));
        if st ~= 0; error('GPT failed (sigma_z %d nm, chunk %d): %s', sz_nm, ic, msg); end
        tg = toc(tc);
        gdf = fullfile(gdir, 'chunk_result.gdf');
        E = zeros(3, numel(lam)); I = E; tr_time = 0;
        for b = 1:2                                               % 2 reads of 1000 particles (~0.7 GB each)
            loc = (b-1)*1000 + (1:1000); tb = tic;
            tr = read_gpt_particles(gdf, loc, true); tr_time = tr_time + toc(tb);
            zp = tr.z(:,1); xp = D0(loc,1) + D0(loc,4)./D0(loc,6).*(zp - D0(loc,3)); yp = D0(loc,2) + D0(loc,5)./D0(loc,6).*(zp - D0(loc,3));
            d = max(abs([tr.x(:,1) - xp; tr.y(:,1) - yp]));
            if d > 1e-9; error('ID / initial distribution mismatch (%.2e m)', d); end
            for k = 1:numel(loc)
                traj = particle_traj(tr, k);
                inmap = abs(traj.r(3,:) - zund) < zmap;
                R.xmax = max(R.xmax, max(abs(traj.r(1:2, inmap)), [], 'all'));
                s = emission_time_spectrum(field_region(traj), obs, lam, 'fB', norm(obs)/c0);
                E = E + s.Ew; I = I + abs(s.Ew).^2;
            end
            clear tr
        end
        delete(gdf);
        R.Ec(ic,:,:) = E; R.Ic(ic,:,:) = I; R.done(ic) = true; R.tchunk(ic) = toc(tc);
        save(out, 'R');
        fprintf('%s sigma_z %2d nm, chunk %2d/%d: GPT %.0f s, read %.0f s, total %.0f s; max |x|,|y| in map %.3f mm  [%s]\n', ...
            tag, sz_nm, ic, nch, tg, tr_time, toc(tc), 1e3*R.xmax, datestr(now, 'HH:MM'));
    end
    if R.xmax > 0.45e-3; warning('a particle reached %.3f mm, close to / beyond the 0.5 mm map', 1e3*R.xmax); end
end
end

function S = standard_sample(res_dir, N)
% standard normal numbers (cut at 4 sigma), the same for every sigma_z
f = fullfile(res_dir, 'standard_sample.mat');
if exist(f, 'file'); S = load(f); return; end
rng(20261104, 'twister'); X = zeros(N, 3); for k = 1:3; X(:, k) = cut_randn(N, 4); end
rng(20261105, 'twister'); U = zeros(N, 3); for k = 1:3; U(:, k) = cut_randn(N, 4); end
E0 = 7e9; me_g = 9.10953e-31; qe_g = 1.6021892e-19; c0 = 299792458;     % GPT constants
G0 = 1 + E0*qe_g/(me_g*c0^2); sx = 20e-6; en = 30e-6; sig_d = 0.02;
S = struct('X', X, 'U', U, 'G0', G0, 'sx', sx, 'sth', en/(G0*sx), 'sig_d', sig_d, 'nmacro', 100e-12/1.602176634e-19/N);
save(f, '-struct', 'S');
fprintf('standard sample: N %d, std of the normals %s, nmacro %.6e\n', N, mat2str(std([X U]), 4), S.nmacro);
end

function D0 = write_chunk(gdir, S, ids, sz, zw)
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

function v = cut_randn(N, ncut)
v = randn(N, 1); out = abs(v) > ncut;
while any(out); v(out) = randn(nnz(out), 1); out = abs(v) > ncut; end
end

function traj = field_region(traj)
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
