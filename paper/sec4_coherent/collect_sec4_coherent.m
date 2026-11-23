%% collect_sec4_coherent
% Collects the multi-electron (2000 electrons, sigma_z = 10 nm, no energy
% spread / emittance) results for the undulator and the dipole into ./data,
% one file per element with all models side by side, and copies the
% original figures into ./figures_original. Nothing is recomputed.
%
% Sources (read only):
%   matlab/radiation_multi_electron/und15p_N2000_onaxis/results/compare_form_factor.mat, compare_translation.mat
%   matlab/radiation_multi_electron/dipole_N2000_onaxis/results/compare_three.mat, check_shot_noise.mat
% Run: matlab -batch "collect_sec4_coherent"

this_dir = fileparts(mfilename('fullpath'));
root     = fullfile(this_dir, '..', '..');
data_dir = fullfile(this_dir, 'data');             if ~exist(data_dir, 'dir'); mkdir(data_dir); end
fig_dir  = fullfile(this_dir, 'figures_original'); if ~exist(fig_dir, 'dir');  mkdir(fig_dir);  end
me_dir   = fullfile(root, 'matlab', 'radiation_multi_electron');

nw_of = @(lam) round(0.01 / mean(diff(log(lam))));          % +-1 % in lambda

%% 1. Undulator
a = load(fullfile(me_dir, 'und15p_N2000_onaxis', 'results', 'compare_form_factor.mat'));  cmp = a.cmp;
b = load(fullfile(me_dir, 'und15p_N2000_onaxis', 'results', 'compare_translation.mat'));  tm  = b.tm;
lam = cmp.lambda_nm(:); omega = cmp.omega(:);
assert(max(abs(tm.lambda_nm(:) - lam) ./ lam) < 1e-14, 'undulator lambda grids differ');
chk = max(abs(tm.S_multi(:) - cmp.S_multi(:))) / max(cmp.S_multi);
assert(chk < 1e-12, 'undulator S_multi differs between the two files (%.2e)', chk);
N = cmp.N;
U = [lam, omega, cmp.S_multi(:), cmp.S1(:), cmp.S_analytic(:), cmp.S_discrete(:), tm.S_translation(:), ...
     cmp.incoherent_sum(:), tm.incoherent_translation(:), cmp.F2_analytic(:), N * cmp.F2_analytic(:), ...
     cmp.err_env_analytic(:), cmp.err_env_discrete(:), tm.err_env_translation(:), cmp.err_env_incoherent(:)];
write_txt(fullfile(data_dir, 'undulator_N2000_onaxis.txt'), U, { ...
    '15-period undulator, 240 MeV, N = 2000 electrons (1 electron per macroparticle), sigma_x = sigma_y = 20 um, sigma_z = 10 nm (sample 9.874 nm),', ...
    'no energy spread, no divergence, no space charge. Observer on axis (0,0,11) m. d2W/domega dOmega [J/(rad/s)/sr] per sr of the 11 m sphere.', ...
    'Models: analytic = S1*[N+N(N-1)|F|^2] (Gaussian, sigma_z = 10 nm); discrete = S1*|sum_j exp(-i w dt_j)|^2 (same sample, phase only);', ...
    'translation = reference trajectory shifted by (x0,y0,0) and in time by -z0/(beta c), exact observation geometry per electron.', ...
    'err_env_* = |MAGNETO - model| / local envelope (max within +-1% in lambda); err_env_incoherent compares sum_j|E_j|^2 with N*S1.', ...
    'Columns: 1 lambda_nm  2 omega_rad_s  3 MAGNETO_multi  4 S1_single_electron  5 analytic_Gaussian  6 discrete_form_factor  7 translation_model', ...
    '  8 incoherent_sum_MAGNETO  9 incoherent_sum_translation  10 |F|^2_analytic  11 N|F|^2  12 err_env_analytic  13 err_env_discrete  14 err_env_translation  15 err_env_incoherent'});

%% 2. Dipole
c = load(fullfile(me_dir, 'dipole_N2000_onaxis', 'results', 'compare_three.mat'));    c3 = c.c3;
d = load(fullfile(me_dir, 'dipole_N2000_onaxis', 'results', 'check_shot_noise.mat')); sn = d.sn;
assert(max(abs(c3.lambda_nm(:) - lam) ./ lam) < 1e-14, 'dipole lambda grid differs from the undulator grid');
assert(max(abs(sn.lambda_nm(:) - lam) ./ lam) < 1e-14, 'shot-noise lambda grid differs');
env = movmax(c3.S_multi(:), 2 * nw_of(lam) + 1);
err_env_incoh = abs(c3.incoherent_multi(:) - c3.N * c3.S1(:)) ./ movmax(c3.N * c3.S1(:), 2 * nw_of(lam) + 1);
Dp = [lam, omega, c3.S_multi(:), c3.S1(:), c3.S_analytic(:), sn.S_discrete(:), c3.S_translation(:), ...
      c3.incoherent_multi(:), c3.incoherent_translation(:), c3.F2(:), c3.N * c3.F2(:), ...
      c3.err_env_analytic(:), sn.err_env_discrete(:), c3.err_env_translation(:), err_env_incoh, ...
      c3.err_env_analytic_sample_sigma(:), c3.S_multi(:) ./ c3.S_analytic(:) - 1];
write_txt(fullfile(data_dir, 'dipole_N2000_onaxis.txt'), Dp, { ...
    'Short dipole 0.1 m / 0.65 T (Enge fringe, gap 1 cm), 240 MeV, N = 2000 electrons, same bunch as the undulator case, tout 3.75 um.', ...
    'Observer on axis (0,0,11) m. d2W/domega dOmega [J/(rad/s)/sr] per sr of the 11 m sphere. Models as in undulator_N2000_onaxis.txt.', ...
    sprintf('err_env_* = |MAGNETO - model| / local envelope (+-1%% in lambda); col 16 uses the sample sigma_z = %.4f nm in the analytic model.', c3.sigma_z_sample * 1e9), ...
    'Columns: 1 lambda_nm  2 omega_rad_s  3 MAGNETO_multi  4 S1_single_electron  5 analytic_Gaussian  6 discrete_form_factor  7 translation_model', ...
    '  8 incoherent_sum_MAGNETO  9 incoherent_sum_translation  10 |F|^2  11 N|F|^2  12 err_env_analytic  13 err_env_discrete  14 err_env_translation', ...
    '  15 err_env_incoherent (sum_j|E_j|^2 vs N*S1)  16 err_env_analytic_sample_sigma  17 MAGNETO/analytic - 1'});

SN = [sn.lambda_ens_nm(:), sn.ens_rel_std(:), sn.sample_rel_dev(:)];
write_txt(fullfile(data_dir, 'dipole_shot_noise.txt'), SN, { ...
    sprintf('Shot noise of a finite Gaussian sample: %d independent samples, N = 2000, truncated at +-4 sigma.', sn.n_realisations), ...
    'ens_rel_std = relative std of |sum_j exp(-i w dt_j)|^2 / [N + N(N-1)|F|^2] over the samples; sample_rel_dev = deviation of the GPT sample.', ...
    'Columns: 1 lambda_nm  2 ens_rel_std  3 sample_rel_dev'});

save(fullfile(data_dir, 'sec4_coherent.mat'), 'U', 'Dp', 'SN', 'N');

%% 3. Key numbers
fid = fopen(fullfile(data_dir, 'summary_numbers.txt'), 'w');
bands = [15 25; 25 40; 40 70; 70 150];
fprintf(fid, '%% max (median) of |MAGNETO - model| / local envelope per band\n');
fprintf(fid, '%% Undulator\n%% band_nm   N|F|^2_range   analytic   discrete   translation   incoherent\n');
for k = 1:4
    j = lam >= bands(k, 1) & lam <= bands(k, 2);
    fprintf(fid, '%3d-%3d  %.2g-%.4g  %s  %s  %s  %s\n', bands(k, :), min(U(j, 11)), max(U(j, 11)), ...
        mm(U(j, 12)), mm(U(j, 13)), mm(U(j, 14)), mm(U(j, 15)));
end
fprintf(fid, '%% Dipole\n%% band_nm   N|F|^2_range   analytic   discrete   translation   incoherent   analytic(sample sigma_z)\n');
for k = 1:4
    j = lam >= bands(k, 1) & lam <= bands(k, 2);
    fprintf(fid, '%3d-%3d  %.2g-%.4g  %s  %s  %s  %s  %s\n', bands(k, :), min(Dp(j, 11)), max(Dp(j, 11)), ...
        mm(Dp(j, 12)), mm(Dp(j, 13)), mm(Dp(j, 14)), mm(Dp(j, 15)), mm(Dp(j, 16)));
end
[~, i150] = min(abs(lam - 150));
fprintf(fid, 'Dipole at 150 nm: S_multi/(N S1) = %.1f, S_multi/(N^2 S1) = %.4f, |F|^2 = %.4f\n', ...
    Dp(i150, 3) / (N * Dp(i150, 4)), Dp(i150, 3) / (N^2 * Dp(i150, 4)), Dp(i150, 10));
for h = [101.8 33.94 20.37]
    [~, ih] = min(abs(lam - h));
    w = max(1, ih - 20):min(numel(lam), ih + 20);          % peak of the harmonic
    [~, ip] = max(U(w, 4)); ip = w(ip);
    fprintf(fid, 'Undulator harmonic near %.2f nm (peak %.2f nm): N|F|^2 = %.3g, S_multi/(N S1) = %.4g, S_multi/(N^2 S1) = %.4g\n', ...
        h, lam(ip), U(ip, 11), U(ip, 3) / (N * U(ip, 4)), U(ip, 3) / (N^2 * U(ip, 4)));
end
fclose(fid);
type(fullfile(data_dir, 'summary_numbers.txt'));

%% 4. Original figures
copy_pngs(fullfile(me_dir, 'und15p_N2000_onaxis', 'figures'), fullfile(fig_dir, 'undulator'));
copy_pngs(fullfile(me_dir, 'dipole_N2000_onaxis', 'figures'), fullfile(fig_dir, 'dipole'));

function s = mm(x)
    s = sprintf('%.2e(%.1e)', max(x), median(x));
end

function write_txt(fn, M, hdr)
    fid = fopen(fn, 'w');
    for k = 1:numel(hdr); fprintf(fid, '%% %s\n', hdr{k}); end
    fmt = [repmat('%.17e ', 1, size(M, 2) - 1) '%.17e\n'];
    fprintf(fid, fmt, M');
    fclose(fid);
end

function copy_pngs(src, dst)
    if ~exist(dst, 'dir'); mkdir(dst); end
    f = dir(fullfile(src, '*.png'));
    for k = 1:numel(f); copyfile(fullfile(src, f(k).name), fullfile(dst, f(k).name)); end
end
