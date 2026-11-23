function run_large_sigz(chunks)
%% run_large_sigz
% User (2026-11-05): SLAC says sigma_z may be 10-1000 nm; two detectors 50-280 nm and 300-1000 nm. Scan sigma_z from
% 10 nm in 10 nm steps until the sensitivity is negligible, at the chosen point x = 0, y = -2.0 mm, z = 5.9 m.
% Same beam line, beam, GPT run (sigma_z = 10 nm), MAGNETO split and dz = 120 um as run_combined.m; one observer,
% lambda 50-1000 nm (600 log points).
% Stored per chunk (no sigma_z dependence, so any sigma_z can be evaluated afterwards):
%   E0(chunk, part, comp, omega) = sum_j E_j shifted to sigma_z = 0 (dt_j = -sz_gpt X3_j / v_zj); part 1 und, 2 dip
%   I (chunk, part, comp, omega) = sum_j |E_j|^2, part 1 und, 2 dip, 3 und + dip
%   Er(chunk, sz, part, comp, omega) = sum_j E_j exp(i w dt_j) for sz = 10:10:100 nm (one macro-particle realisation,
%       for validation of the expected-value formula in analyze_large_sigz.m)
% Resumable: results/large_sigz.mat saved after every chunk; progress in results/large_sigz_progress.log.

this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
% (read_gpt_particles and particle_traj are in magneto/)                                      % read_gpt_particles, particle_traj
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole');
res_dir = fullfile(this_dir, 'results');
maxNumCompThreads(8);

c0 = 299792458; dz_um = 120; zsplit = 0.63; sz_gpt = 10e-9;
obs = [0 -2e-3 5.9];
lam = logspace(log10(50), log10(1000), 600)*1e-9; om = 2*pi*c0./lam; nw = numel(lam);
szr = (10:10:100)*1e-9; nsr = numel(szr);
Ntot = 5e4; Nc = 2000; nch = Ntot/Nc;
S = load(fullfile(this_dir, '..', 'E7GeV_gap6und_sigz_sweep', 'results', 'standard_sample.mat'));
S.X(:, 1:2) = -S.X(:, 1:2); S.U(:, 1:2) = -S.U(:, 1:2);           % rotation by 180 deg about z (as run_combined)

out = fullfile(res_dir, 'large_sigz.mat'); logf = fullfile(res_dir, 'large_sigz_progress.log');
if exist(out, 'file'); R = load(out); R = R.R;
else
    R = struct('sz_gpt', sz_gpt, 'obs', obs, 'lambda', lam, 'omega', om, 'dz_um', dz_um, 'sz_real', szr, ...
        'E0', zeros(nch, 2, 3, nw), 'I', zeros(nch, 3, 3, nw), 'Er', zeros(nch, nsr, 2, 3, nw), 'done', false(1, nch), ...
        'nmacro', S.nmacro, 'tchunk', nan(1, nch), ...
        'note', ['E0 = sum_j E_j at sigma_z = 0 (time shift), I = sum_j |E_j|^2, Er = realisation at sz_real; ' ...
                 'part 1 undulator (z < 0.63 m), 2 dipole (3 = und + dip for I); comp x y z; fields carry q*nmacro']);
end
for ic = chunks
    if R.done(ic); continue; end
    tc = tic; ids = (ic-1)*Nc + (1:Nc);
    D0 = write_chunk(gdir, S, ids, sz_gpt, 0.4);
    [st, msg] = system(sprintf('cmd /c ""%s" %d"', fullfile(gdir, 'run.bat'), dz_um));
    if st ~= 0; error('GPT failed (chunk %d): %s', ic, msg); end
    tg = toc(tc);
    gdf = fullfile(gdir, 'chunk_result.gdf');
    E0 = zeros(2, 3, nw); I = zeros(3, 3, nw); Er = zeros(nsr, 2, 3, nw);
    X3 = S.X(ids, 3); vz = c0*D0(:,6)./D0(:,7);
    for b = 1:2
        loc = (b-1)*1000 + (1:1000);
        tr = read_gpt_particles(gdf, loc, true);
        zp = tr.z(:,1); xp = D0(loc,1) + D0(loc,4)./D0(loc,6).*(zp - D0(loc,3)); yp = D0(loc,2) + D0(loc,5)./D0(loc,6).*(zp - D0(loc,3));
        if max(abs([tr.x(:,1) - xp; tr.y(:,1) - yp])) > 1e-9; error('ID / initial distribution mismatch'); end
        for k = 1:numel(loc)
            j = loc(k); traj = particle_traj(tr, k);
            gap = traj.r(3,:) > 0.6245 & traj.r(3,:) < 0.6390;
            if any(traj.fB(:, gap) ~= 0, 'all'); error('field between the magnets is not zero'); end
            Eu = emission_time_spectrum(sub(traj, traj.r(3,:) < zsplit), obs, lam, 'fB', norm(obs)/c0).Ew;
            Ed = emission_time_spectrum(sub(traj, traj.r(3,:) >= zsplit), obs, lam, 'fB', norm(obs)/c0).Ew;
            I(1,:,:) = squeeze(I(1,:,:)) + abs(Eu).^2;
            I(2,:,:) = squeeze(I(2,:,:)) + abs(Ed).^2;
            I(3,:,:) = squeeze(I(3,:,:)) + abs(Eu + Ed).^2;
            p0 = exp(1i*(-sz_gpt*X3(j)/vz(j))*om);                      % shift to sigma_z = 0
            E0(1,:,:) = squeeze(E0(1,:,:)) + Eu.*p0;
            E0(2,:,:) = squeeze(E0(2,:,:)) + Ed.*p0;
            ph = exp(1i*((szr(:) - sz_gpt)*X3(j)/vz(j))*om);            % nsr x nw
            for s = 1:nsr
                Er(s,1,:,:) = squeeze(Er(s,1,:,:)) + Eu.*ph(s,:);
                Er(s,2,:,:) = squeeze(Er(s,2,:,:)) + Ed.*ph(s,:);
            end
        end
        clear tr
    end
    delete(gdf);
    R.E0(ic,:,:,:) = E0; R.I(ic,:,:,:) = I; R.Er(ic,:,:,:,:) = Er; R.done(ic) = true; R.tchunk(ic) = toc(tc);
    save(out, 'R');
    fid = fopen(logf, 'a');
    fprintf(fid, 'chunk %2d/%d: GPT %.0f s, total %.0f s  [%s]\n', ic, nch, tg, toc(tc), datestr(now, 'HH:MM'));
    fclose(fid);
end
end

function D0 = write_chunk(gdir, S, ids, sz, zw)
% identical to run_combined.m / run_sweep.m
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
k = find(m); kf = k(any(t.fB(:, k) ~= 0, 1)); k = max(k(1), kf(1)-1) : min(k(end), kf(end)+1);
for f = {'t', 'G'};             t.(f{1}) = t.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; t.(f{1}) = t.(f{1})(:, k); end
end
