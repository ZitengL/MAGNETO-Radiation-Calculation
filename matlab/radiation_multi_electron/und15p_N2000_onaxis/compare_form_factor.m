%% compare_form_factor
% Benchmark of the MAGNETO multi-particle result (run_multi_onaxis.m) against the
% single-electron spectrum combined with a form factor, on axis at z = 11 m.
% Without space charge, energy spread or divergence every electron follows the
% reference trajectory shifted by its initial offset, so its field is the
% single-electron field with a phase factor:
%   E~_j(w) = E~_1(w) exp(-i w dt_j),   dt_j = -z0_j/(beta c) + (x0_j^2 + y0_j^2)/(2 L c)
% (z0 ahead -> earlier arrival; the second term is the longer path of an off-axis
% electron to the on-axis observer, L = 11 m - 1 m = 10 m from the undulator centre).
%  (1) discrete form factor (exact test, same sampled particles):
%        S_disc(w) = S_1(w) |sum_j exp(-i w dt_j)|^2
%  (2) analytic expectation for a Gaussian bunch:
%        S_an(w) = S_1(w) [N + N(N-1) |F(w)|^2],  F = exp(-(w sigma_z / c)^2 / 2)
%      (the transverse factor 1/(1+(k sigma_perp^2/L)^2) per plane is > 0.9997 here and is neglected)
%  (3) incoherent check: sum_j |E~_j|^2 = N S_1
% Outputs: results/compare_form_factor.(mat|txt), figures/compare_form_factor.(png|fig),
%          figures/compare_form_factor_ratio.(png|fig)

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
c0 = 299792458; eps0 = 8.85418781e-12;

M = load(fullfile(res_dir, 'multi_onaxis.mat')); M = M.multi;
lam = M.lambda_nm; w = M.omega; N = M.params.N; obs = M.params.obs_point_m; t_ref = M.params.t_ref_s;
pref = eps0 * c0 * norm(obs)^2 / pi;

%% single-electron reference (same field map, same output grid)
traj = read_gpt_trajectory(fullfile(root_dir, 'gpt_sims', 'single_electron', 'E240MeV_und15p_paraxial', 'result.gdf'));
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
s1 = emission_time_spectrum(traj, obs, lam*1e-9, 'fB', t_ref);
S1 = pref * sum(abs(s1.Ew).^2, 1);
h = harmonic_wavelengths(traj, [1 3 5]);
beta = sqrt(1 - 1/traj.G(1)^2);

%% form factors
D0 = readmatrix(fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E240MeV_und15p_N2000', 'initial_distribution.txt'), 'NumHeaderLines', 1);
x0 = D0(:,1); y0 = D0(:,2); z0 = D0(:,3);
L = norm(obs) - 1.0;
dt = -z0/(beta*c0) + (x0.^2 + y0.^2)/(2*L*c0);
B = zeros(1, numel(w));
for j0 = 1:500:numel(w)                                   % |sum_j exp(-i w dt_j)|^2 in blocks
    jj = j0:min(j0+499, numel(w));
    B(jj) = abs(sum(exp(-1i * dt * w(jj)), 1)).^2;
end
sig_z = 10e-9;
F2 = exp(-(w*sig_z/c0).^2);
S_disc = S1 .* B;
S_an   = S1 .* (N + N*(N-1)*F2);
S_mult = M.d2W_dwdOmega.total;
S_inc  = M.incoherent_sum.total;

%% numbers
% Pointwise ratios blow up at the spectral nulls between the sidelobes (S1 ~ 0), so the
% differences are measured relative to the local envelope: moving maximum over +-1% in lambda.
nw = round(0.01 / mean(diff(log(lam))));
env = @(S) movmax(S, 2*nw + 1);
e_disc = abs(S_mult - S_disc) ./ env(S_disc);
e_an   = abs(S_mult - S_an)   ./ env(S_an);
e_inc  = abs(S_inc - N*S1)    ./ env(N*S1);
r_disc = S_mult ./ S_disc - 1;                            % pointwise (for the record)
fprintf('pointwise multi/disc - 1 at the S1 nulls reaches %.1e; below: differences relative to the local envelope (+-1%% in lambda)\n', max(abs(r_disc)));
fprintf('incoherent check |sum_j|E_j|^2 - N S1| / envelope: max %.3e\n', max(e_inc));
for band = [15 25; 25 40; 40 70; 70 150]'
    sel = lam >= band(1) & lam < band(2);
    fprintf('  %3d-%3d nm: |multi - disc|/env max %.2e (median %.1e); |multi - analytic|/env max %.2e; incoherent check %.1e; N|F|^2 %.2e .. %.2e\n', band, ...
        max(e_disc(sel)), median(e_disc(sel)), max(e_an(sel)), max(e_inc(sel)), N*min(F2(sel)), N*max(F2(sel)));
end
lam_h = h.lambda_beta*1e9;
fprintf('harmonics (beta method) n = 1, 3, 5: %s nm; N|F|^2 there: %s\n', mat2str(lam_h, 5), mat2str(N*exp(-(2*pi*sig_z./(lam_h*1e-9)).^2), 3));
for i = 1:3
    [~, j] = min(abs(lam - lam_h(i)));
    win = abs(lam - lam(j)) < 0.03*lam(j); [~, jj] = max(S1 .* win);
    fprintf('  n = %d peak (%.2f nm): S_multi / (N S1) = %.1f, S_multi / (N^2 S1) = %.3f, |F|^2 = %.3e\n', ...
        2*i-1, lam(jj), S_mult(jj)/(N*S1(jj)), S_mult(jj)/(N^2*S1(jj)), F2(jj));
end

%% save
cmp.lambda_nm = lam; cmp.omega = w; cmp.S1 = S1; cmp.S_multi = S_mult; cmp.S_discrete = S_disc; cmp.S_analytic = S_an;
cmp.incoherent_sum = S_inc; cmp.bunching_discrete = B; cmp.err_env_discrete = e_disc; cmp.err_env_analytic = e_an; cmp.err_env_incoherent = e_inc; cmp.F2_analytic = F2; cmp.harmonics_nm = lam_h; cmp.N = N;
cmp.note = 'all spectra d2W/domega dOmega [J/(rad/s)/sr], on axis at (0,0,11) m; bunching_discrete = |sum_j exp(-i w dt_j)|^2';
save(fullfile(res_dir, 'compare_form_factor.mat'), 'cmp');
fid = fopen(fullfile(res_dir, 'compare_form_factor.txt'), 'wt');
fprintf(fid, '%% On axis at (0,0,11) m, d2W/domega dOmega [J/(rad/s)/sr]; N = %d electrons, sigma_z = 10 nm\n', N);
fprintf(fid, '%% Columns: 1 lambda [nm]  2 omega [rad/s]  3 MAGNETO multi-particle (coherent sum)  4 single electron S1  5 S1*|sum_j exp(-i w dt_j)|^2 (discrete)  6 S1*[N+N(N-1)|F|^2] (analytic)  7 sum_j|E_j|^2 (incoherent sum)\n');
fprintf(fid, '%-14.8f %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e\n', [lam; w; S_mult; S1; S_disc; S_an; S_inc]);
fclose(fid);

%% figures
xt = [15 20 30 40 50 70 100 150];
f = figure('Position', [100 100 950 600]);
loglog(lam, S_mult, 'k-', 'LineWidth', 1.6); hold on;
loglog(lam, S_disc, 'r--', 'LineWidth', 1.1);
loglog(lam, S_an, 'b:', 'LineWidth', 1.3);
loglog(lam, N*S1, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
loglog(lam, N^2*S1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
for i = 1:3; xline(lam_h(i), 'm:', sprintf('n = %d', 2*i-1), 'LabelVerticalAlignment', 'bottom'); end
grid on; xlim([15 150]); set(gca, 'XTick', xt);
legend('MAGNETO, 2000 electrons (coherent sum)', 'single electron \times discrete form factor', ...
       'single electron \times [N + N(N-1)|F|^2] (Gaussian)', 'N \times single (incoherent)', 'N^2 \times single (fully coherent)', ...
       'Location', 'southeast');
xlabel('\lambda [nm]'); ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]');
title('Coherent radiation benchmark: 15-period undulator, 240 MeV, \sigma_z = 10 nm, on axis at 11 m');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'compare_form_factor.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'compare_form_factor.fig'));

f = figure('Position', [100 100 950 600]);
loglog(lam, e_disc, 'r-', lam, e_an, 'b-', lam, e_inc, 'k-'); grid on; xlim([15 150]); set(gca, 'XTick', xt);
legend('|multi - single \times discrete form factor|', '|multi - single \times [N+N(N-1)|F|^2]|', '|\Sigma|E_j|^2 - N S_1|', 'Location', 'best');
xlabel('\lambda [nm]'); ylabel('difference / local envelope (\pm1% in \lambda)');
title('MAGNETO multi-particle vs single electron \times form factor');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'compare_form_factor_ratio.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'compare_form_factor_ratio.fig'));
