function verify_timeshift_gpt(chunks)
%% verify_timeshift_gpt
% User (2026-11-05): check the time-shift method against direct GPT runs at 3 sigma_z (user's choice: 100, 300, 500 nm;
% 1e4 particles = chunks 1-5 of the standard sample).
% For each chunk, GPT is run at sigma_z = 10 nm (the reference used in run_large_sigz.m) and at 100, 300, 500 nm,
% with the same random numbers (only z = sigma_z * X3 changes). For every particle j:
%   direct  : E_j(sz) from the GPT run at sz
%   shifted : E_j(10 nm) * exp(i w (sz - 10 nm) X3_j / v_zj)     (the method used in run_large_sigz.m)
% Same observer, wavelengths, dz, undulator/dipole split as run_large_sigz.m.
% Stored per chunk (resumable, results/verify_timeshift_gpt.mat):
%   Sd(chunk, s, part, comp, w), Ss(...)   : sum_j of direct / shifted fields (s = 100, 300, 500 nm; part 1 und, 2 dip)
%   Id, Is (chunk, s, comp, w)             : sum_j |E_und + E_dip|^2, direct / shifted
%   E0(chunk, r, part, comp, w)            : fields shifted to sz = 0 from each GPT run r = 10, 100, 300, 500 nm
%   I0(chunk, r, comp, w)                  : sum_j |E_und + E_dip|^2 of each GPT run
%   dpp(particle, s, band)                 : per-particle || E_direct - E_shifted || / || E_direct || (x comp, und + dip,
%                                            L2 over the frequencies of band 1 D1 50-280 nm, 2 D2 300-1000 nm)
% Progress: results/verify_timeshift_progress.log

this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
% (read_gpt_particles and particle_traj are in magneto/)
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole');
res_dir = fullfile(this_dir, 'results');
maxNumCompThreads(8);

c0 = 299792458; dz_um = 120; zsplit = 0.63;
obs = [0 -2e-3 5.9];
lam = logspace(log10(50), log10(1000), 600)*1e-9; om = 2*pi*c0./lam; nw = numel(lam);   % as run_large_sigz.m
szr = [10 100 300 500]*1e-9; nr = numel(szr); ns = nr - 1;                              % run 1 = reference
bands = [50e-9 280e-9; 300e-9 1000e-9]; mb = {lam >= bands(1,1) & lam <= bands(1,2), lam >= bands(2,1) & lam <= bands(2,2)};
Nc = 2000; nchmax = 5;
S = load(fullfile(this_dir, '..', 'E7GeV_gap6und_sigz_sweep', 'results', 'standard_sample.mat'));
S.X(:, 1:2) = -S.X(:, 1:2); S.U(:, 1:2) = -S.U(:, 1:2);           % rotation by 180 deg about z (as run_large_sigz)

out = fullfile(res_dir, 'verify_timeshift_gpt.mat'); logf = fullfile(res_dir, 'verify_timeshift_progress.log');
if exist(out, 'file'); V = load(out); V = V.V;
else
    V = struct('obs', obs, 'lambda', lam, 'omega', om, 'dz_um', dz_um, 'sz_runs', szr, 'bands', bands, ...
        'Sd', zeros(nchmax, ns, 2, 3, nw), 'Ss', zeros(nchmax, ns, 2, 3, nw), 'Id', zeros(nchmax, ns, 3, nw), ...
        'Is', zeros(nchmax, ns, 3, nw), 'E0', zeros(nchmax, nr, 2, 3, nw), 'I0', zeros(nchmax, nr, 3, nw), ...
        'dpp', nan(nchmax*Nc, ns, 2), 'done', false(1, nchmax), 'nmacro', S.nmacro, 'tchunk', nan(1, nchmax));
end
for ic = chunks
    if V.done(ic); continue; end
    tc = tic; ids = (ic-1)*Nc + (1:Nc); X3 = S.X(ids, 3);
    Eref = zeros(Nc, 2, 3, nw);                                    % per-particle fields of the 10 nm run
    for r = 1:nr
        D0 = write_chunk(gdir, S, ids, szr(r), 0.4);
        vz = c0*D0(:,6)./D0(:,7);
        gdf = fullfile(gdir, 'chunk_result.gdf'); if exist(gdf, 'file'); delete(gdf); end
        [st, msg] = system(sprintf('cmd /c ""%s" %d"', fullfile(gdir, 'run.bat'), dz_um));
        if st ~= 0 || ~exist(gdf, 'file'); error('GPT failed (chunk %d, sz %g): %s', ic, szr(r), msg); end
        for b = 1:2
            loc = (b-1)*1000 + (1:1000);
            tr = read_gpt_particles(gdf, loc, true);
            zp = tr.z(:,1); xp = D0(loc,1) + D0(loc,4)./D0(loc,6).*(zp - D0(loc,3)); yp = D0(loc,2) + D0(loc,5)./D0(loc,6).*(zp - D0(loc,3));
            if max(abs([tr.x(:,1) - xp; tr.y(:,1) - yp])) > 1e-9; error('ID / initial distribution mismatch'); end
            for k = 1:numel(loc)
                j = loc(k); traj = particle_traj(tr, k);
                Eu = emission_time_spectrum(sub(traj, traj.r(3,:) < zsplit), obs, lam, 'fB', norm(obs)/c0).Ew;
                Ed = emission_time_spectrum(sub(traj, traj.r(3,:) >= zsplit), obs, lam, 'fB', norm(obs)/c0).Ew;
                p0 = exp(1i*(-szr(r)*X3(j)/vz(j))*om);              % shift this run to sigma_z = 0
                V.E0(ic,r,1,:,:) = squeeze(V.E0(ic,r,1,:,:)) + Eu.*p0;
                V.E0(ic,r,2,:,:) = squeeze(V.E0(ic,r,2,:,:)) + Ed.*p0;
                V.I0(ic,r,:,:) = squeeze(V.I0(ic,r,:,:)) + abs(Eu + Ed).^2;
                if r == 1
                    Eref(j,1,:,:) = Eu; Eref(j,2,:,:) = Ed;
                else
                    s = r - 1; ph = exp(1i*((szr(r) - szr(1))*X3(j)/vz(j))*om);
                    Esu = squeeze(Eref(j,1,:,:)).*ph; Esd = squeeze(Eref(j,2,:,:)).*ph;
                    V.Sd(ic,s,1,:,:) = squeeze(V.Sd(ic,s,1,:,:)) + Eu;  V.Sd(ic,s,2,:,:) = squeeze(V.Sd(ic,s,2,:,:)) + Ed;
                    V.Ss(ic,s,1,:,:) = squeeze(V.Ss(ic,s,1,:,:)) + Esu; V.Ss(ic,s,2,:,:) = squeeze(V.Ss(ic,s,2,:,:)) + Esd;
                    V.Id(ic,s,:,:) = squeeze(V.Id(ic,s,:,:)) + abs(Eu + Ed).^2;
                    V.Is(ic,s,:,:) = squeeze(V.Is(ic,s,:,:)) + abs(Esu + Esd).^2;
                    ed = Eu(1,:) + Ed(1,:); es = Esu(1,:) + Esd(1,:);
                    for bb = 1:2
                        V.dpp(ids(j), s, bb) = norm(ed(mb{bb}) - es(mb{bb}))/norm(ed(mb{bb}));
                    end
                end
            end
            clear tr
        end
        delete(gdf);
    end
    V.done(ic) = true; V.tchunk(ic) = toc(tc);
    save(out, 'V');
    fid = fopen(logf, 'a'); fprintf(fid, 'chunk %d/%d: %.0f s  [%s]\n', ic, nchmax, toc(tc), datestr(now, 'HH:MM')); fclose(fid);
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
