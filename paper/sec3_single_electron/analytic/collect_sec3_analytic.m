%% collect_sec3_analytic
% Collects the single-electron analytic comparisons (undulator: Wiedemann
% 4th ed. (26.46)/(26.51); dipole: Jackson (14.79)) into ./data, copies the
% original figures into ./figures_original. Nothing is recomputed.
%
% Sources (read only):
%   matlab/undulator_analytic_comparison/undulator240/results/compare_undulator.mat
%   matlab/undulator_analytic_comparison/undulator240/diagnostic_textbook_trajectory/results/diag_textbook.mat
%   matlab/dipole_analytic_comparison/dipole240_tangent/results/compare_dipole_tangent.mat (+ txt)
%   matlab/dipole_analytic_comparison/dipole240_onaxis/results/compare_dipole_onaxis.txt
% Run: matlab -batch "collect_sec3_analytic"

this_dir = fileparts(mfilename('fullpath'));
root     = fullfile(this_dir, '..', '..', '..');
data_dir = fullfile(this_dir, 'data');             if ~exist(data_dir, 'dir'); mkdir(data_dir); end
fig_dir  = fullfile(this_dir, 'figures_original'); if ~exist(fig_dir, 'dir');  mkdir(fig_dir);  end
und_dir  = fullfile(root, 'matlab', 'undulator_analytic_comparison', 'undulator240');
dip_dir  = fullfile(root, 'matlab', 'dipole_analytic_comparison');

%% 1. Undulator: analytic, MAGNETO ideal field, MAGNETO real map, MAGNETO on the textbook trajectory
a = load(fullfile(und_dir, 'results', 'compare_undulator.mat'));          out = a.out;
b = load(fullfile(und_dir, 'diagnostic_textbook_trajectory', 'results', 'diag_textbook.mat')); D = b.D;
lam = out.lambda_nm(:); omega = out.omega(:);
assert(max(abs(D.lambda_nm(:) - lam) ./ lam) < 1e-14, 'lambda grids differ');

obs_names = {'11m_theta0', '11m_theta0p5', '1km_theta0', '1km_theta0p5'};
fid_s = fopen(fullfile(data_dir, 'summary_numbers.txt'), 'w');
fprintf(fid_s, '%% Undulator vs analytic (Wiedemann 26.46), per sr seen from the undulator centre\n');
und = struct();
for k = 1:numel(obs_names)
    c = out.(['c_' obs_names{k}]);
    t = D.dz375.(['c_' obs_names{k}]);
    chk = max(abs(t.S_analytic(:) - c.S_analytic_2646(:))) / max(c.S_analytic_2646);
    assert(chk < 1e-10, 'analytic spectra of the two sources differ (%s, %.2e)', obs_names{k}, chk);
    M = [lam, omega, c.S_analytic_2646(:), c.S_analytic_2651(:), c.S_MAGNETO_ideal(:), c.S_MAGNETO_real(:), ...
         t.S_textbook_traj(:), c.err_env_ideal(:), c.err_env_real(:), t.err_env(:)];
    und.(['c_' obs_names{k}]) = struct('M', M, 'obs', c.obs, 'theta_u', c.theta_u, 'R_src', c.R_src, ...
        'peaks', c.peaks, 'peaks_columns', c.peaks_columns);
    fn = fullfile(data_dir, sprintf('undulator_analytic_%s.txt', obs_names{k}));
    fid = fopen(fn, 'w');
    fprintf(fid, '%% 240 MeV single electron, 15-period undulator (K = %.6f, lambda_u = %.6f mm, N = 15 in the analytic formula)\n', ...
        out.params.K, out.params.lambda_u * 1e3);
    fprintf(fid, '%% Observer %s: obs = (%.6g, %.6g, %.6g) m, theta_u = %.4f mrad (seen from the undulator centre), R_src = %.4f m\n', ...
        obs_names{k}, c.obs, c.theta_u * 1e3, c.R_src);
    fprintf(fid, '%% d2W/domega dOmega [J/(rad/s)/sr] per sr seen from the undulator centre, total polarisation\n');
    fprintf(fid, '%% textbook trajectory = hard-edged cosine field with exactly 15 periods (the trajectory assumed by the formula), MAGNETO, sampling 3.75 um\n');
    fprintf(fid, '%% err_env = (MAGNETO - analytic 26.46) / local envelope of the analytic spectrum (max within +-1%% in lambda)\n');
    fprintf(fid, '%% Columns: 1 lambda_nm  2 omega_rad_s  3 analytic_26.46  4 analytic_26.51(near peaks only)  5 MAGNETO_ideal_unduplan  6 MAGNETO_real_map  7 MAGNETO_textbook_trajectory  8 err_env_ideal  9 err_env_real  10 err_env_textbook\n');
    fprintf(fid, '%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n', M');
    fclose(fid);
    % textbook trajectory vs analytic, where the analytic spectrum > 1% of its maximum
    sel = c.S_analytic_2646(:) > 0.01 * max(c.S_analytic_2646);
    fprintf(fid_s, '%s: max |MAGNETO(textbook trajectory) - analytic| / local envelope, where S > 1%% of max: %.2e\n', ...
        obs_names{k}, max(abs(t.err_env(sel))));
    fprintf(fid_s, '   peaks table, columns: %s\n', c.peaks_columns);
    fprintf(fid_s, '   %s\n', strjoin(compose('%.4g', c.peaks(1, :)), ' '));
    for r = 2:size(c.peaks, 1); fprintf(fid_s, '   %s\n', strjoin(compose('%.4g', c.peaks(r, :)), ' ')); end
end

%% 2. Dipole on the arc-midpoint tangent: analytic, hard-edge MAGNETO, Enge-fringe MAGNETO
fprintf(fid_s, '%%\n%% Dipole on the arc-midpoint tangent (per sr seen from the tangent point)\n');
dip = struct();
c = load(fullfile(dip_dir, 'dipole240_tangent', 'results', 'compare_dipole_tangent.mat'));  ct = c.out;
assert(max(abs(ct.lambda_nm(:) - lam) ./ lam) < 1e-14, 'dipole tangent lambda grid differs');
for tag = {'11m', '1km'}
    o = ct.(['o' tag{1}]);
    r_hard = o.ratio_hard_analytic_m1(:);  r_enge = o.ratio_enge_analytic_m1(:);
    M = [lam, omega, o.S_analytic(:), o.S_MAGNETO_hardedge(:), o.S_MAGNETO_enge(:), r_hard, r_enge];
    dip.(['tangent_' tag{1}]) = M;
    fid = fopen(fullfile(data_dir, ['dipole_tangent_analytic_' tag{1} '.txt']), 'w');
    fprintf(fid, '%% 240 MeV single electron, dipole 0.1 m / 0.65 T (bend %.4f mrad, rho = %.6f m), observer on the arc-midpoint tangent at (%.6g, %.6g, %.6g) m\n', ...
        ct.bend_rad * 1e3, ct.rho, o.obs);
    fprintf(fid, '%% d2W/domega dOmega [J/(rad/s)/sr] per sr seen from the tangent point (R_src = %.4f m), total polarisation; full precision from compare_dipole_tangent.mat\n', o.R_src);
    fprintf(fid, '%% hard edge: b1 = 1e5 1/m; Enge: b1 = 200 1/m (gap 1 cm), the dipole used everywhere else in the paper\n');
    fprintf(fid, '%% Columns: 1 lambda_nm  2 omega_rad_s  3 analytic Jackson (14.79) theta = 0  4 MAGNETO hard-edge dipole  5 MAGNETO Enge-fringe dipole  6 hard/analytic-1  7 Enge/analytic-1\n');
    fprintf(fid, '%.17e %.17e %.17e %.17e %.17e %.17e %.17e\n', M');
    fclose(fid);
    fprintf(fid_s, '%s: hard/analytic-1 in [%+.2e, %+.2e]; Enge/analytic-1 in [%+.2e, %+.2e]\n', ...
        tag{1}, min(r_hard), max(r_hard), min(r_enge), max(r_enge));
end

%% 3. Dipole on axis vs the same formula (formula not applicable there; kept for the explanation)
copyfile(fullfile(dip_dir, 'dipole240_onaxis', 'results', 'compare_dipole_onaxis.txt'), ...
         fullfile(data_dir, 'dipole_onaxis_vs_analytic.txt'));
A = readmatrix(fullfile(data_dir, 'dipole_onaxis_vs_analytic.txt'), 'FileType', 'text', 'CommentStyle', '%');
dip.onaxis = A;
fprintf(fid_s, '%%\n%% Dipole on axis: MAGNETO / analytic at 15, 50, 150 nm = %.2e, %.3f, %.3f (edge radiation, formula not applicable)\n', ...
    interp1(A(:, 1), A(:, 5), 15), interp1(A(:, 1), A(:, 5), 50), interp1(A(:, 1), A(:, 5), 150));
fclose(fid_s);

save(fullfile(data_dir, 'sec3_analytic.mat'), 'und', 'dip', 'lam', 'omega');

%% 4. Original figures (png) for reference
copy_pngs(fullfile(und_dir, 'figures'), fullfile(fig_dir, 'undulator'));
copy_pngs(fullfile(und_dir, 'diagnostic_textbook_trajectory', 'figures'), fullfile(fig_dir, 'undulator_textbook_trajectory'));
copy_pngs(fullfile(dip_dir, 'dipole240_tangent', 'figures'), fullfile(fig_dir, 'dipole_tangent'));
copy_pngs(fullfile(dip_dir, 'dipole240_onaxis', 'figures'), fullfile(fig_dir, 'dipole_onaxis'));
type(fullfile(data_dir, 'summary_numbers.txt'));

function copy_pngs(src, dst)
    if ~exist(dst, 'dir'); mkdir(dst); end
    f = dir(fullfile(src, '*.png'));
    for k = 1:numel(f); copyfile(fullfile(src, f(k).name), fullfile(dst, f(k).name)); end
end
