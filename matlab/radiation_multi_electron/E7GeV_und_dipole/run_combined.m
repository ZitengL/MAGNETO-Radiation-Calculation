function run_combined(chunks, sz_gpt_nm, szl_nm, tag)
%% run_combined
% SLAC gap-6 undulator (turned over) + dipole (bent up, 0.13 T, Enge gap 35 mm) tracked together, 7 GeV kinetic,
% 100 pC as 5e4 macro-particles, no space charge. User's decisions (2026-11-05):
%   track undulator and dipole together; observers in the dipole bend plane x = 0 on the side opposite to the bend,
%   y = -1, -1.5, -2, -2.5, -3 mm at z = 5.9 m (x polariser removes the dipole light there by symmetry);
%   sigma_z = 10:5:50 nm by time shift this time (a full direct check later); energy in 50-300 nm.
% Beam: the same standard normal numbers as E7GeV_gap6und_sigz_sweep (results/standard_sample.mat), rotated by
%   180 deg about z (x, y, x', y' change sign) so that, with the undulator turned over, the undulator part is the
%   exact rotated image of the earlier sigma_z sweep (validation of the combined run and of the time shift).
% Time shift: without space charge a particle started dz further along z follows the same path (transverse shift
%   x' dz ~ 1e-11 m, neglected) dt = dz / v_z earlier, so E_j(sz) = E_j(sz_gpt) exp(i w dt_j),
%   dz_j = (sz - sz_gpt) X3_j (emission_time_spectrum uses exp(-i w (t_obs - t_ref))).
% MAGNETO: emission_time_spectrum, t_ref = |obs|/c, fields carry q*nmacro; the integral is split at z = 0.63 m
%   (field-free between 0.6238 and 0.6394 m) into the undulator and the dipole part.
% Stored per chunk: Ec(chunk, sz, obs, part, comp, omega) = sum_j E_j (part 1 undulator, 2 dipole) and
%   Ic(chunk, obs, part, comp, omega) = sum_j |E_j|^2 (part 1 und, 2 dip, 3 und + dip; independent of sz).
% Resumable: results/<tag>combined.mat is saved after every chunk; progress in results/<tag>progress.log.
%
% Syntax: run_combined(1:25, 10, 10:5:50, '')      main run (GPT at sigma_z = 10 nm, time shift to 15..50 nm)
%         run_combined(1:3, 50, 50, 'direct50_')   direct GPT check at sigma_z = 50 nm

this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
% (read_gpt_particles and particle_traj are in magneto/)                                      % read_gpt_particles, particle_traj
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole');
res_dir = fullfile(this_dir, 'results'); if ~exist(res_dir, 'dir'); mkdir(res_dir); end
maxNumCompThreads(12);

c0 = 299792458; dz_um = 120; zsplit = 0.63; zund = 0.4; zmap = 0.22381;
ys = [-1 -1.5 -2 -2.5 -3]*1e-3; obs = [zeros(numel(ys),1) ys(:) 5.9*ones(numel(ys),1)]; no = size(obs, 1);
lam = logspace(log10(50), log10(300), 300)*1e-9; om = 2*pi*c0./lam; nw = numel(lam);
Ntot = 5e4; Nc = 2000; nch = Ntot/Nc; nsz = numel(szl_nm);
S = load(fullfile(this_dir, '..', 'E7GeV_gap6und_sigz_sweep', 'results', 'standard_sample.mat'));
S.X(:, 1:2) = -S.X(:, 1:2); S.U(:, 1:2) = -S.U(:, 1:2);           % rotation by 180 deg about z

out = fullfile(res_dir, [tag 'combined.mat']); logf = fullfile(res_dir, [tag 'progress.log']);
if exist(out, 'file'); R = load(out); R = R.R;
else
    R = struct('sz_gpt', sz_gpt_nm*1e-9, 'sigma_z', szl_nm*1e-9, 'obs', obs, 'lambda', lam, 'omega', om, 'dz_um', dz_um, ...
        'Ec', zeros(nch, nsz, no, 2, 3, nw), 'Ic', zeros(nch, no, 3, 3, nw), 'done', false(1, nch), ...
        'xmax_und', 0, 'ymin_end', inf, 'ymax_end', -inf, 'nmacro', S.nmacro, 'tchunk', nan(1, nch), ...
        'note', ['Ec(chunk,sz,obs,part,comp,omega) = sum_j E_j, part 1 undulator (z<0.63 m), 2 dipole; Ic(chunk,obs,part,comp,omega) = sum_j |E_j|^2, ' ...
                 'part 1 und, 2 dip, 3 und+dip; MAGNETO = eps0 c |obs|^2/pi |sum_chunks Ec|^2; beam = standard sample rotated by 180 deg about z']);
end
for ic = chunks
    if R.done(ic); continue; end
    tc = tic; ids = (ic-1)*Nc + (1:Nc);
    D0 = write_chunk(gdir, S, ids, sz_gpt_nm*1e-9, 0.4);
    [st, msg] = system(sprintf('cmd /c ""%s" %d"', fullfile(gdir, 'run.bat'), dz_um));
    if st ~= 0; error('GPT failed (chunk %d): %s', ic, msg); end
    tg = toc(tc);
    gdf = fullfile(gdir, 'chunk_result.gdf');
    Ec = zeros(nsz, no, 2, 3, nw); Ic = zeros(no, 3, 3, nw); tr_time = 0;
    X3 = S.X(ids, 3); vz = c0*D0(:,6)./D0(:,7);                     % initial v_z
    for b = 1:2                                                     % 2 reads of 1000 particles
        loc = (b-1)*1000 + (1:1000); tb = tic;
        tr = read_gpt_particles(gdf, loc, true); tr_time = tr_time + toc(tb);
        zp = tr.z(:,1); xp = D0(loc,1) + D0(loc,4)./D0(loc,6).*(zp - D0(loc,3)); yp = D0(loc,2) + D0(loc,5)./D0(loc,6).*(zp - D0(loc,3));
        d = max(abs([tr.x(:,1) - xp; tr.y(:,1) - yp]));
        if d > 1e-9; error('ID / initial distribution mismatch (%.2e m)', d); end
        for k = 1:numel(loc)
            j = loc(k); traj = particle_traj(tr, k);
            inmap = abs(traj.r(3,:) - zund) < zmap;
            R.xmax_und = max(R.xmax_und, max(abs(traj.r(1:2, inmap)), [], 'all'));
            yend = traj.r(2,end) + traj.b(2,end)/traj.b(3,end)*(5.9 - traj.r(3,end));
            R.ymin_end = min(R.ymin_end, yend); R.ymax_end = max(R.ymax_end, yend);
            gap = traj.r(3,:) > 0.6245 & traj.r(3,:) < 0.6390;            % undulator map ends 0.62381, dipole fringe starts 0.63943
            if any(traj.fB(:, gap) ~= 0, 'all'); error('field between the magnets is not zero'); end
            seg = {sub(traj, traj.r(3,:) < zsplit), sub(traj, traj.r(3,:) >= zsplit)};
            dt = (szl_nm(:)*1e-9 - sz_gpt_nm*1e-9)*X3(j)/vz(j);                % nsz x 1
            ph = exp(1i*dt*om);                                                % nsz x nw
            for o = 1:no
                Eu = emission_time_spectrum(seg{1}, obs(o,:), lam, 'fB', norm(obs(o,:))/c0).Ew;
                Ed = emission_time_spectrum(seg{2}, obs(o,:), lam, 'fB', norm(obs(o,:))/c0).Ew;
                Ic(o,1,:,:) = squeeze(Ic(o,1,:,:)) + abs(Eu).^2;
                Ic(o,2,:,:) = squeeze(Ic(o,2,:,:)) + abs(Ed).^2;
                Ic(o,3,:,:) = squeeze(Ic(o,3,:,:)) + abs(Eu + Ed).^2;
                for s = 1:nsz
                    Ec(s,o,1,:,:) = squeeze(Ec(s,o,1,:,:)) + Eu.*ph(s,:);
                    Ec(s,o,2,:,:) = squeeze(Ec(s,o,2,:,:)) + Ed.*ph(s,:);
                end
            end
        end
        clear tr
    end
    delete(gdf);
    R.Ec(ic,:,:,:,:,:) = Ec; R.Ic(ic,:,:,:,:) = Ic; R.done(ic) = true; R.tchunk(ic) = toc(tc);
    save(out, 'R');
    fid = fopen(logf, 'a');
    fprintf(fid, '%s chunk %2d/%d: GPT %.0f s, read %.0f s, total %.0f s; max |x|,|y| in undulator map %.3f mm; y at 5.9 m %.2f..%.2f mm  [%s]\n', ...
        tag, ic, nch, tg, tr_time, toc(tc), 1e3*R.xmax_und, 1e3*R.ymin_end, 1e3*R.ymax_end, datestr(now, 'HH:MM'));
    fclose(fid);
end
end

function D0 = write_chunk(gdir, S, ids, sz, zw)
% same formulas as E7GeV_gap6und_sigz_sweep/run_sweep.m (straight-line drift from the start to the waist at zw)
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
