%% run_und240_spectra
% Single-electron spectra of the 15-period undulator at 240 MeV (near-axis paraxial map),
% observer on the mirror sphere R = 11 m in the wiggle plane: theta_x = 0, 0.1, 0.2, 0.3, 0.4,
% 0.5 mrad (+x side), D = 11 m (sin theta_x, 0, cos theta_x). 15-150 nm, 5000 log points
% (same grid as the multi-electron benchmark). Energy per unit angular frequency.
% Trajectory: gpt_sims/single_electron/E240MeV_und15p_paraxial/result.gdf (tout 60 um),
% convergence check with the 30 and 15 um runs (tout_convergence/).
% Harmonic markers: lambda_n(theta) = lambda_1(0) [1 + gamma^2 theta_u^2/(1 + K^2/2)] / n,
% lambda_1(0) and K from the <beta_z> method (harmonic_wavelengths.m); theta_u is the angle seen
% from the undulator centre (z = 1 m), theta_u = atan(11 sin(theta)/(11 cos(theta) - 1)) ~ 1.1 theta.
% On axis only odd harmonics are marked (even ones vanish), off axis n = 1..5.
%
% Outputs: results/und240_spectra.mat, results/und240_spectra_energy.txt, results/und240_spectra_photons.txt
%          figures/und240_theta<k>_{logy,liny}.(png|fig) (k = 0..5 -> 0.0..0.5 mrad),
%          figures/und240_overlay.(png|fig), figures/und240_tout_convergence.(png|fig)

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(this_dir));
addpath(fullfile(root_dir, 'magneto'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
if ~exist(res_dir, 'dir'); mkdir(res_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
c0 = 299792458; hbar = 1.05457182e-34;

run_dir = fullfile(root_dir, 'gpt_sims', 'single_electron', 'E240MeV_und15p_paraxial');
files = {fullfile(run_dir, 'result.gdf'), fullfile(run_dir, 'tout_convergence', 'result_dz30um.gdf'), ...
         fullfile(run_dir, 'tout_convergence', 'result_dz15um.gdf')};
dz_um = [60 30 15];
R = 11; zu = 1.0;
theta = (0:0.1:0.5) * 1e-3;
lambda_nm = logspace(log10(15), log10(150), 5000); lambda = lambda_nm * 1e-9;

%% trajectories and harmonics
T = cell(1, 3);
for f = 1:3; T{f} = field_region(read_gpt_trajectory(files{f})); end
h = harmonic_wavelengths(read_gpt_trajectory(files{1}), 1:5);
gam = T{1}.G(1); K = h.K_from_beta; lam1 = h.lambda_beta(1) * 1e9;
theta_u = atan(R*sin(theta) ./ (R*cos(theta) - zu));
fprintf('gamma %.4f, lambda_u %.4f mm, K (<beta_z>) %.4f, lambda_1(0) = %.3f nm\n', gam, h.lambda_u*1e3, K, lam1);

%% spectra
nt = numel(theta); nl = numel(lambda);
W = zeros(nt, nl); Wc = zeros(3, nt, nl); Wpol = zeros(nt, 3, nl);
for i = 1:nt
    D = R * [sin(theta(i)), 0, cos(theta(i))];
    for f = 1:3
        s = emission_time_spectrum(T{f}, D, lambda, 'fB');
        Wc(f, i, :) = s.d2W_dwdOmega.total;
        if f == 1
            W(i, :) = s.d2W_dwdOmega.total;
            Wpol(i, :, :) = [s.d2W_dwdOmega.x; s.d2W_dwdOmega.y; s.d2W_dwdOmega.z];
        end
    end
end
Np = W ./ (hbar * 2*pi*c0 ./ lambda);                  % photons per unit angular frequency

%% convergence (relative to the local envelope, +-1% in lambda)
nw = round(0.01 / mean(diff(log(lambda_nm)))); env = @(S) movmax(S, 2*nw + 1);
conv_err = zeros(2, nt);
for i = 1:nt
    s15 = squeeze(Wc(3, i, :))';
    conv_err(1, i) = max(abs(squeeze(Wc(1, i, :))' - s15) ./ env(s15));
    conv_err(2, i) = max(abs(squeeze(Wc(2, i, :))' - s15) ./ env(s15));
    fprintf('theta_x = %.1f mrad (theta_u = %.3f mrad): max |S60 - S15|/env = %.2e, |S30 - S15|/env = %.2e; peak %.3e J/(rad/s)/sr at %.2f nm\n', ...
        theta(i)*1e3, theta_u(i)*1e3, conv_err(1, i), conv_err(2, i), max(W(i, :)), lambda_nm(W(i, :) == max(W(i, :))));
end

%% harmonic positions
lam_h = lam1 * (1 + gam^2 * theta_u'.^2 / (1 + K^2/2)) ./ (1:5);   % nt x 5

%% save
sp.lambda_nm = lambda_nm; sp.omega = 2*pi*c0 ./ lambda; sp.theta_x_rad = theta; sp.theta_u_rad = theta_u;
sp.d2W_dwdOmega = W; sp.d2W_dwdOmega_xyz = Wpol; sp.d2N_dwdOmega = Np; sp.harmonics_nm = lam_h;
sp.conv_dz_um = dz_um; sp.conv_d2W_dwdOmega = Wc; sp.conv_err_env = conv_err;
sp.note = 'rows = theta_x; d2W/domega dOmega [J/(rad/s)/sr] per sr of the 11 m sphere centred at the origin; xyz = polarisation components';
save(fullfile(res_dir, 'und240_spectra.mat'), 'sp');
for q = 1:2
    if q == 1; fn = 'und240_spectra_energy.txt'; Q = W; u = 'd2W/domega dOmega [J/(rad/s)/sr]';
    else;      fn = 'und240_spectra_photons.txt'; Q = Np; u = 'd2N/domega dOmega [photons/(rad/s)/sr]'; end
    fid = fopen(fullfile(res_dir, fn), 'wt');
    fprintf(fid, '%% 15-period undulator, single electron, 240 MeV, observer on the R = 11 m sphere in the wiggle plane, total polarisation\n');
    fprintf(fid, '%% Columns: 1 lambda [nm]  2 omega [rad/s]  3-8 %s at theta_x = 0, 0.1, 0.2, 0.3, 0.4, 0.5 mrad\n', u);
    fprintf(fid, ['%-14.8f %-14.6e' repmat(' %-14.6e', 1, nt) '\n'], [lambda_nm; sp.omega; Q]);
    fclose(fid);
end

%% figures
xt = [15 20 30 40 50 70 100 150];
for i = 1:nt
    for m = 1:2
        f = figure('Position', [100 100 900 500]);
        if m == 1; loglog(lambda_nm, W(i, :), 'k-', 'LineWidth', 1.1); tag = 'logy';
        else;      semilogx(lambda_nm, W(i, :), 'k-', 'LineWidth', 1.1); tag = 'liny'; end
        hold on; grid on; xlim([15 150]); set(gca, 'XTick', xt);
        if theta(i) == 0; nlist = [1 3 5]; else; nlist = 1:5; end
        for n = nlist
            xline(lam_h(i, n), 'r--', sprintf('n = %d: %.2f nm', n, lam_h(i, n)), 'LineWidth', 1, ...
                'LabelOrientation', 'aligned', 'LabelVerticalAlignment', 'middle', 'FontSize', 8);
        end
        xlabel('\lambda [nm]'); ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]');
        title({sprintf('15-period undulator, 240 MeV, single electron, \\theta_x = %.1f mrad (R = 11 m)', theta(i)*1e3), ...
               'red dashed: harmonics \lambda_n(\theta) = \lambda_1(0)[1 + \gamma^2\theta_u^2/(1+K^2/2)]/n'});
        set(findall(f, 'Type', 'axes'), 'Toolbar', []);
        exportgraphics(f, fullfile(fig_dir, sprintf('und240_theta%d_%s.png', i-1, tag)), 'Resolution', 200);
        savefig(f, fullfile(fig_dir, sprintf('und240_theta%d_%s.fig', i-1, tag)));
        close(f);
    end
end
f = figure('Position', [100 100 950 600]);
cols = lines(nt);
for i = 1:nt; loglog(lambda_nm, W(i, :), '-', 'Color', cols(i, :), 'LineWidth', 1.1); hold on; end
grid on; xlim([15 150]); set(gca, 'XTick', xt);
legend(arrayfun(@(t) sprintf('\\theta_x = %.1f mrad', t*1e3), theta, 'UniformOutput', false), 'Location', 'southeast');
xlabel('\lambda [nm]'); ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]');
title('15-period undulator, 240 MeV, single electron: spectra at six angles in the wiggle plane');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'und240_overlay.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'und240_overlay.fig'));
f = figure('Position', [100 100 950 500]);
for i = 1:nt
    s15 = squeeze(Wc(3, i, :))';
    loglog(lambda_nm, abs(squeeze(Wc(1, i, :))' - s15) ./ env(s15), '-', 'Color', cols(i, :)); hold on;
end
grid on; xlim([15 150]); set(gca, 'XTick', xt);
legend(arrayfun(@(t) sprintf('\\theta_x = %.1f mrad', t*1e3), theta, 'UniformOutput', false), 'Location', 'best');
xlabel('\lambda [nm]'); ylabel('|S_{60\mum} - S_{15\mum}| / local envelope');
title('tout step convergence (60 vs 15 \mum)');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'und240_tout_convergence.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'und240_tout_convergence.fig'));

function traj = field_region(traj)
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
