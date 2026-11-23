%% check_shot_noise
% Check of the MAGNETO - analytic Gaussian difference in the incoherent / transition region
% (15-40 nm). Not part of the three published results.
%  (1) discrete form factor of the same 2000 electrons: S_disc = S1 |sum_j exp(-i w dt_j)|^2,
%      dt_j = -z0/(beta c) + (x0^2 + y0^2)/(2 L c), L = 10 m (as in ../und15p_N2000_onaxis)
%      -> if MAGNETO = S_disc, the difference from the analytic curve is the shot noise of this sample.
%  (2) expected scatter: 200 independent Gaussian samples (N = 2000, sigma_z = 10 nm, sigma_x,y = 20 um,
%      truncated at +-4 sigma by rejection), relative standard deviation of
%      |sum exp(-i w dt_j)|^2 / [N + N(N-1)|F|^2] at every 10th wavelength.
% Output: results/check_shot_noise.mat, figures/check_shot_noise.(png|fig)

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; sig_z = 10e-9; sig_t = 20e-6; L = 10;

C = load(fullfile(res_dir, 'compare_three.mat')); C = C.c3;
lam = C.lambda_nm; w = C.omega; N = C.N; S1 = C.S1;
M = load(fullfile(res_dir, 'multi_onaxis.mat')); gam = 470.66442480; beta = sqrt(1 - 1/gam^2);
D0 = readmatrix(fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E240MeV_dipole_N2000', 'initial_distribution.txt'), 'NumHeaderLines', 1);
bunch = @(x, y, z, ww) abs(sum(exp(-1i * (-z/(beta*c0) + (x.^2 + y.^2)/(2*L*c0)) * ww), 1)).^2;

%% (1) discrete form factor of this sample
B = zeros(1, numel(w));
for j0 = 1:500:numel(w); jj = j0:min(j0+499, numel(w)); B(jj) = bunch(D0(:,1), D0(:,2), D0(:,3), w(jj)); end
S_disc = S1 .* B;
nw = round(0.01 / mean(diff(log(lam)))); env = @(S) movmax(S, 2*nw + 1);
e_disc = abs(C.S_multi - S_disc) ./ env(S_disc);
fprintf('|MAGNETO - S1 x discrete form factor| / env: max %.2e over 15-150 nm\n', max(e_disc));

%% (2) ensemble scatter
rng(1); nr = 200; iw = 1:10:numel(w); R = zeros(nr, numel(iw));
G = N + N*(N-1)*exp(-(w(iw)*sig_z/c0).^2);
for r = 1:nr
    X = zeros(N, 3); n = 0;
    while n < N
        u = randn(N, 3) .* [sig_t sig_t sig_z]; u = u(all(abs(u) <= 4*[sig_t sig_t sig_z], 2), :);
        k = min(size(u, 1), N - n); X(n+1:n+k, :) = u(1:k, :); n = n + k;
    end
    R(r, :) = bunch(X(:,1), X(:,2), X(:,3), w(iw)) ./ G;
end
sd = std(R, 0, 1);
dev = C.S_multi(iw) ./ C.S_analytic(iw) - 1;                % this sample, pointwise relative to the analytic
for band = [15 25; 25 40; 40 70; 70 150]'
    sel = lam(iw) >= band(1) & lam(iw) <= band(2);
    fprintf('  %3d-%3d nm: expected rel. std of the shot noise %.2f..%.2f; this sample |S_multi/S_an - 1| max %.3f, max |dev|/std %.2f\n', ...
        band, min(sd(sel)), max(sd(sel)), max(abs(dev(sel))), max(abs(dev(sel)) ./ sd(sel)));
end

sn.lambda_nm = lam; sn.S_discrete = S_disc; sn.err_env_discrete = e_disc; sn.lambda_ens_nm = lam(iw);
sn.ens_rel_std = sd; sn.sample_rel_dev = dev; sn.n_realisations = nr;
sn.note = 'S_discrete = S1 |sum_j exp(-i w dt_j)|^2 of the GPT sample; ens_rel_std = std over random Gaussian samples of B/[N+N(N-1)|F|^2]';
save(fullfile(res_dir, 'check_shot_noise.mat'), 'sn');

f = figure('Position', [100 100 950 500]);
semilogx(lam(iw), dev, 'k-', 'LineWidth', 1.3); hold on;
semilogx(lam(iw), sd, 'b--', lam(iw), -sd, 'b--', lam(iw), 2*sd, 'b:', lam(iw), -2*sd, 'b:');
grid on; xlim([15 150]); set(gca, 'XTick', [15 20 30 40 50 70 100 150]);
legend('this sample: S_{MAGNETO} / S_{analytic} - 1', '\pm1 std of the shot noise (200 Gaussian samples)', '', '\pm2 std', 'Location', 'northeast');
xlabel('\lambda [nm]'); ylabel('relative deviation');
title('Dipole, 2000 electrons: deviation from the analytic Gaussian form factor vs expected shot noise');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'check_shot_noise.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'check_shot_noise.fig'));
