%% compare_three
% Coherent-radiation benchmark for the short dipole (240 MeV, 0.1 m, 0.65 T), on axis at 11 m:
%   (1) MAGNETO multi-particle: coherent sum over the 2000 GPT trajectories (run_multi_onaxis.m)
%   (2) analytic Gaussian form factor: S1 [N + N(N-1)|F|^2], F = exp(-(w sigma_z/c)^2/2), sigma_z = 10 nm
%   (3) strict translation model: the single-electron reference trajectory shifted by
%       (x0, y0, 0) and in time by -z0/(beta c) for every electron, summed coherently
%       (same as compare_translation.m of the undulator benchmark)
% S1 = single-electron spectrum from gpt_sims/single_electron/E240MeV_dipole_L10cm (same dipole,
% same output grid, tout 3.75 um). Accuracy: lambda_acc = smallest wavelength above which the
% on-axis single-electron spectrum from the 3.75 um trajectory agrees with the 1.875 um one to 1%
% (pointwise); with 3.75 um this is the whole 15-150 nm range (tout_sweep.m). (With the earlier
% 60 um step it was only ~65 nm and above, see results_dz60um/.) The txt files and figures
% contain only lambda >= lambda_acc.
% Outputs: results/dipole_multi_MAGNETO.txt, results/dipole_analytic_formfactor.txt,
%          results/dipole_translation.txt, results/compare_three.mat,
%          figures/compare_three.(png|fig), figures/compare_three_difference.(png|fig)

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
c0 = 299792458; eps0 = 8.85418781e-12; sig_z = 10e-9;

M = load(fullfile(res_dir, 'multi_onaxis.mat')); M = M.multi;
lam = M.lambda_nm; w = M.omega; N = M.params.N; obs = M.params.obs_point_m; t_ref = M.params.t_ref_s;
pref = eps0 * c0 * norm(obs)^2 / pi;
sdir = fullfile(root_dir, 'gpt_sims', 'single_electron', 'E240MeV_dipole_L10cm');

%% single-electron reference and its accuracy range
ref = field_region(read_gpt_trajectory(fullfile(sdir, 'result.gdf')));
s1 = emission_time_spectrum(ref, obs, lam*1e-9, 'fB', t_ref);
S1 = pref * sum(abs(s1.Ew).^2, 1);
r15 = field_region(read_gpt_trajectory(fullfile(sdir, 'tout_convergence', 'result_dz1p875um.gdf')));
s15 = emission_time_spectrum(r15, obs, lam*1e-9, 'fB', t_ref);
S15 = pref * sum(abs(s15.Ew).^2, 1);                     % finer reference (1.875 um)
i_acc = max([0, find(abs(S1 - S15) > 0.01 * S15, 1, 'last')]) + 1;
lambda_acc = lam(i_acc); ok = i_acc:numel(lam);
fprintf('single electron, 3.75 vs 1.875 um: pointwise within 1%% for lambda >= %.3f nm (max |S/S_fine - 1| = %.2e)\n', ...
    lambda_acc, max(abs(S1./S15 - 1)));

%% (2) analytic Gaussian form factor
F2 = exp(-(w*sig_z/c0).^2);
S_an = S1 .* (N + N*(N-1)*F2);

%% (3) translation model
beta0 = sqrt(1 - 1/ref.G(1)^2);
D0 = readmatrix(fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E240MeV_dipole_N2000', 'initial_distribution.txt'), 'NumHeaderLines', 1);
Esum = zeros(3, numel(w)); incoh = zeros(1, numel(w));
tic;
for j = 1:N
    tr = ref;
    tr.t = ref.t - D0(j, 3) / (beta0 * c0);
    tr.r = ref.r + [D0(j, 1); D0(j, 2); 0];
    s = emission_time_spectrum(tr, obs, lam*1e-9, 'fB', t_ref);
    Esum = Esum + s.Ew; incoh = incoh + sum(abs(s.Ew).^2, 1);
end
fprintf('translation model: %.1f min\n', toc/60);
S_tr = pref * sum(abs(Esum).^2, 1);

%% (1) MAGNETO multi
S_mul = M.d2W_dwdOmega.total;

%% differences (relative to the local envelope, +-1% in lambda; pointwise ratios blow up at nulls)
nw = round(0.01 / mean(diff(log(lam)))); env = @(S) movmax(S, 2*nw + 1);
e_an = abs(S_mul - S_an) ./ env(S_an);
e_tr = abs(S_mul - S_tr) ./ env(S_tr);
e_inc = abs(M.incoherent_sum.total - N*S1) ./ env(N*S1);
fprintf('converged range %.2f-150 nm, N|F|^2 = %.0f .. %.0f\n', lambda_acc, N*min(F2(ok)), N*max(F2(ok)));
sz_s = std(D0(:, 3)); S_an_s = S1 .* (N + N*(N-1)*exp(-(w*sz_s/c0).^2));   % check only: sample sigma_z
e_an_s = abs(S_mul - S_an_s) ./ env(S_an_s);
for band = [max(lambda_acc, 15) 25; 25 40; 40 70; 70 150]'
    sel = lam >= band(1) & lam <= band(2) & lam >= lambda_acc;
    if ~any(sel); continue; end
    fprintf('  %6.2f-%3.0f nm: N|F|^2 %.2e..%.2e; analytic with sample sigma_z (%.4f nm): max %.2e\n', band, N*min(F2(sel)), N*max(F2(sel)), sz_s*1e9, max(e_an_s(sel)));
    fprintf('  %6.2f-%3.0f nm: |multi - analytic|/env max %.2e (median %.1e); |multi - translation|/env max %.2e (median %.1e); |sum|E_j|^2 - N S1|/env max %.1e\n', ...
        band, max(e_an(sel)), median(e_an(sel)), max(e_tr(sel)), median(e_tr(sel)), max(e_inc(sel)));
end
fprintf('at 150 nm: S_multi / (N S1) = %.1f, S_multi / (N^2 S1) = %.3f, |F|^2 = %.4f\n', S_mul(end)/(N*S1(end)), S_mul(end)/(N^2*S1(end)), F2(end));

%% save
c3.lambda_nm = lam; c3.omega = w; c3.S1 = S1; c3.S1_dz1p875um = S15; c3.S_multi = S_mul; c3.S_analytic = S_an; c3.S_translation = S_tr;
c3.incoherent_multi = M.incoherent_sum.total; c3.incoherent_translation = pref*incoh; c3.F2 = F2; c3.lambda_acc_nm = lambda_acc;
c3.err_env_analytic = e_an; c3.err_env_translation = e_tr; c3.N = N; c3.sigma_z_sample = sz_s; c3.err_env_analytic_sample_sigma = e_an_s;
c3.note = 'd2W/domega dOmega [J/(rad/s)/sr], on axis (0,0,11) m; full 15-150 nm grid, converged only for lambda >= lambda_acc_nm';
save(fullfile(res_dir, 'compare_three.mat'), 'c3');
hdr = sprintf('%% Short dipole (240 MeV, 0.1 m, 0.65 T), 2000 electrons, sigma_z = 10 nm, on axis at (0,0,11) m. Converged range only (lambda >= %.3f nm).\n', lambda_acc);
files = {'dipole_multi_MAGNETO.txt', 'dipole_analytic_formfactor.txt', 'dipole_translation.txt'};
desc  = {'MAGNETO multi-particle: coherent sum over the GPT trajectories', ...
         'single electron x analytic Gaussian form factor: S1 [N + N(N-1) exp(-(w sigma_z/c)^2)], sigma_z = 10 nm', ...
         'translation model: single-electron trajectory shifted by (x0,y0,0) and -z0/(beta c), coherent sum'};
data  = {S_mul, S_an, S_tr};
for q = 1:3
    fid = fopen(fullfile(res_dir, files{q}), 'wt');
    fprintf(fid, '%s%% %s\n%% Columns: 1 lambda [nm]  2 omega [rad/s]  3 d2W/domega dOmega [J/(rad/s)/sr]\n', hdr, desc{q});
    fprintf(fid, '%-14.8f %-14.6e %-14.6e\n', [lam(ok); w(ok); data{q}(ok)]);
    fclose(fid);
end

%% figures (converged range)
xt = [15 20 30 40 50 70 100 150]; xt = xt(xt >= lambda_acc); L = lam(ok);
f = figure('Position', [100 100 950 600]);
loglog(L, S_mul(ok), 'k-', 'LineWidth', 1.6); hold on;
loglog(L, S_an(ok), 'b--', 'LineWidth', 1.2); loglog(L, S_tr(ok), 'r:', 'LineWidth', 1.4);
loglog(L, N*S1(ok), 'Color', [0.5 0.5 0.5]); loglog(L, N^2*S1(ok), '--', 'Color', [0.5 0.5 0.5]);
grid on; xlim([lambda_acc 150]); set(gca, 'XTick', xt);
legend('MAGNETO, 2000 electrons (coherent sum)', 'single electron \times [N + N(N-1)|F|^2] (Gaussian)', ...
       'translation model', 'N \times single', 'N^2 \times single', 'Location', 'best');
xlabel('\lambda [nm]'); ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]');
title(sprintf('Coherent benchmark, dipole 0.1 m / 0.65 T, 240 MeV, \\sigma_z = 10 nm, on axis at 11 m (\\lambda \\geq %.1f nm)', lambda_acc));
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'compare_three.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'compare_three.fig'));
f = figure('Position', [100 100 950 500]);
loglog(L, e_an(ok), 'b-', L, e_tr(ok), 'r-', L, e_inc(ok), 'k-'); grid on; xlim([lambda_acc 150]); set(gca, 'XTick', xt);
legend('|MAGNETO - analytic Gaussian|', '|MAGNETO - translation model|', '|\Sigma|E_j|^2 - N S_1|', 'Location', 'best');
xlabel('\lambda [nm]'); ylabel('difference / local envelope (\pm1% in \lambda)');
title('Dipole: MAGNETO multi-particle vs models');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'compare_three_difference.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'compare_three_difference.fig'));

function traj = field_region(traj)
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
