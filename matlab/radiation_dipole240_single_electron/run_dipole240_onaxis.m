%% run_dipole240_onaxis
% Single-electron spectrum of the short dipole (240 MeV, centre z = 1 m, 0.1 m, 0.65 T,
% bend ~81 mrad, lambda_c ~50 nm), observer on axis at (0, 0, 11) m, 15-150 nm, 5000 log
% points (same grid as the multi-electron benchmark). Energy per unit angular frequency.
% Trajectory: gpt_sims/single_electron/E240MeV_dipole_L10cm/result.gdf (tout 3.75 um);
% convergence with the 7.5 and 1.875 um runs. The on-axis observer sees the entrance of the bend;
% the rest of the arc points up to 81 mrad away, and its emission phase needs a fine step at
% short wavelengths: with 60 um the result converged only above ~65 nm (tout_sweep.m), with
% 3.75 um over the whole 15-150 nm. Only lambda >= lambda_acc (3.75 vs 1.875 um within 1%
% pointwise) would be shown in the figures and txt; the .mat keeps the full grid.
% Outputs: results/dipole240_onaxis.mat, results/dipole240_onaxis.txt,
%          figures/dipole240_onaxis_{logy,liny}.(png|fig), figures/dipole240_tout_convergence.(png|fig)

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(this_dir));
addpath(fullfile(root_dir, 'magneto'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
if ~exist(res_dir, 'dir'); mkdir(res_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
c0 = 299792458; hbar = 1.05457182e-34;
run_dir = fullfile(root_dir, 'gpt_sims', 'single_electron', 'E240MeV_dipole_L10cm');
files = {fullfile(run_dir, 'result.gdf'), fullfile(run_dir, 'tout_convergence', 'result_dz7p5um.gdf'), ...
         fullfile(run_dir, 'tout_convergence', 'result_dz1p875um.gdf')};
dz_lab = {'3.75 \mum (main)', '7.5 \mum', '1.875 \mum'};
obs = [0 0 11];
lambda_nm = logspace(log10(15), log10(150), 5000); lambda = lambda_nm * 1e-9;

S = zeros(3, numel(lambda));
for f = 1:3
    tr = field_region(read_gpt_trajectory(files{f}));
    s = emission_time_spectrum(tr, obs, lambda, 'fB');
    S(f, :) = s.d2W_dwdOmega.total;
    if f == 1; s60 = s; end
end
W = S(1, :); Wpol = [s60.d2W_dwdOmega.x; s60.d2W_dwdOmega.y; s60.d2W_dwdOmega.z];
Np = W ./ (hbar * 2*pi*c0 ./ lambda);
i_acc = max([0, find(abs(S(1,:) - S(3,:)) > 0.01 * S(3,:), 1, 'last')]) + 1; lambda_acc = lambda_nm(i_acc); ok = i_acc:numel(lambda);
j75 = max([0, find(abs(S(2,:) - S(3,:)) > 0.01 * S(3,:), 1, 'last')]) + 1;
fprintf('3.75 vs 1.875 um within 1%% (pointwise) for lambda >= %.3f nm; 7.5 vs 1.875 um: >= %.3f nm\n', lambda_acc, lambda_nm(j75));
fprintf('in the converged range: max |S_3.75/S_1.875 - 1| = %.2e; S(lambda_acc)/max = %.3e; max at %.2f nm\n', ...
    max(abs(S(1,ok)./S(3,ok) - 1)), W(i_acc)/max(W(ok)), lambda_nm(ok(W(ok) == max(W(ok)))));
fprintf('polarisation: max x/total %.1e, max z/total %.1e\n', max(Wpol(1,:)./W), max(Wpol(3,:)./W));

d.lambda_nm = lambda_nm; d.omega = 2*pi*c0 ./ lambda; d.d2W_dwdOmega = W; d.d2W_dwdOmega_xyz = Wpol; d.d2N_dwdOmega = Np;
d.conv_dz_um = [3.75 7.5 1.875]; d.conv_d2W_dwdOmega = S; d.lambda_acc_nm = lambda_acc;
d.note = 'on axis (0,0,11) m, d2W/domega dOmega [J/(rad/s)/sr]; full grid, converged only for lambda >= lambda_acc_nm';
save(fullfile(res_dir, 'dipole240_onaxis.mat'), 'd');
fid = fopen(fullfile(res_dir, 'dipole240_onaxis.txt'), 'wt');
fprintf(fid, '%% Short dipole (240 MeV, 0.1 m, 0.65 T), single electron, on axis at (0,0,11) m. Converged range only (lambda >= %.3f nm).\n', lambda_acc);
fprintf(fid, '%% Columns: 1 lambda [nm]  2 omega [rad/s]  3-6 d2W/domega dOmega [J/(rad/s)/sr] total,x,y,z  7 d2N/domega dOmega [photons/(rad/s)/sr] total\n');
fprintf(fid, '%-14.8f %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e\n', [lambda_nm(ok); d.omega(ok); W(ok); Wpol(:, ok); Np(ok)]);
fclose(fid);

xt = [15 20 30 40 50 70 100 150]; xt = xt(xt >= lambda_acc);
for m = 1:2
    f = figure('Position', [100 100 900 500]);
    if m == 1; loglog(lambda_nm(ok), W(ok), 'k-', 'LineWidth', 1.2); tag = 'logy';
    else;      semilogx(lambda_nm(ok), W(ok), 'k-', 'LineWidth', 1.2); tag = 'liny'; end
    grid on; xlim([lambda_acc 150]); set(gca, 'XTick', xt);
    xlabel('\lambda [nm]'); ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]');
    title({'Short dipole (0.1 m, 0.65 T, bend 81 mrad, \lambda_c = 49.6 nm), 240 MeV, single electron, on axis at 11 m', ...
           sprintf('converged range \\lambda \\geq %.1f nm', lambda_acc)});
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(f, fullfile(fig_dir, ['dipole240_onaxis_' tag '.png']), 'Resolution', 200);
    savefig(f, fullfile(fig_dir, ['dipole240_onaxis_' tag '.fig']));
end
f = figure('Position', [100 100 900 550]);
subplot(2,1,1); loglog(lambda_nm, S(1,:), 'k-', lambda_nm, S(2,:), 'b--', lambda_nm, S(3,:), 'r:'); grid on; xlim([15 150]);
legend(dz_lab{:}, 'Location', 'best');
ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]'); title('Dipole 240 MeV, on axis: tout step convergence (full grid; full sweep in tout_sweep.png)', 'Interpreter', 'none');
subplot(2,1,2); loglog(lambda_nm, abs(S(1,:)./S(3,:) - 1), 'k-', lambda_nm, abs(S(2,:)./S(3,:) - 1), 'b--'); grid on; xlim([15 150]);
yline(0.01, 'k:'); legend('|3.75/1.875 - 1|', '|7.5/1.875 - 1|', 'Location', 'best'); xlabel('\lambda [nm]'); ylabel('relative difference');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'dipole240_tout_convergence.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'dipole240_tout_convergence.fig'));

function traj = field_region(traj)
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
