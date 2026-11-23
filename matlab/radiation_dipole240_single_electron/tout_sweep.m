%% tout_sweep
% tout-step convergence of the on-axis dipole spectrum (240 MeV, 0.1 m, 0.65 T) down to
% 1.875 um: dz = 60, 30, 15, 7.5, 3.75, 1.875 um (gpt_sims/single_electron/E240MeV_dipole_L10cm,
% tout_convergence/). For each step: lambda_acc = smallest wavelength above which the result
% agrees with the next finer step and with the finest one to 1% (pointwise), and the max
% pointwise difference over 15-150 nm. Also the cost (field-region samples).
% Output: results/tout_sweep.mat, figures/tout_sweep.(png|fig)

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(this_dir));
addpath(fullfile(root_dir, 'magneto'));
d = fullfile(root_dir, 'gpt_sims', 'single_electron', 'E240MeV_dipole_L10cm');
dz = [60 30 15 7.5 3.75 1.875];
f = {fullfile(d, 'tout_convergence', 'result_dz60um.gdf'), fullfile(d, 'tout_convergence', 'result_dz30um.gdf'), fullfile(d, 'tout_convergence', 'result_dz15um.gdf'), ...
     fullfile(d, 'tout_convergence', 'result_dz7p5um.gdf'), fullfile(d, 'tout_convergence', 'result_dz3p75um.gdf'), ...
     fullfile(d, 'tout_convergence', 'result_dz1p875um.gdf')};
lam = logspace(log10(15), log10(150), 5000);
S = zeros(numel(dz), numel(lam)); ns = zeros(1, numel(dz)); tc = zeros(1, numel(dz));
for i = 1:numel(dz)
    t = read_gpt_trajectory(f{i});
    k = find(any(t.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(t.t), k(end)+1);
    for g = {'t', 'G'}; t.(g{1}) = t.(g{1})(k); end
    for g = {'r', 'b', 'fE', 'fB'}; t.(g{1}) = t.(g{1})(:, k); end
    tic; s = emission_time_spectrum(t, [0 0 11], lam*1e-9, 'fB', 11/299792458); tc(i) = toc;
    S(i, :) = s.d2W_dwdOmega.total; ns(i) = numel(t.t);
end
lacc = @(a, b) lam(min(numel(lam), max([0, find(abs(a - b) > 0.01*b, 1, 'last')]) + 1));
fprintf('dz [um] | samples | time [s] | lambda_acc vs next finer | vs finest | max |S/S_finest - 1|\n');
for i = 1:numel(dz)
    if i < numel(dz); ln = lacc(S(i,:), S(i+1,:)); else; ln = NaN; end
    fprintf('%6.3f  %7d  %6.2f  %8.3f nm  %8.3f nm  %.2e\n', dz(i), ns(i), tc(i), ln, lacc(S(i,:), S(end,:)), max(abs(S(i,:)./S(end,:) - 1)));
end
save(fullfile(this_dir, 'results', 'tout_sweep.mat'), 'dz', 'lam', 'S', 'ns', 'tc');
fg = figure('Position', [100 100 950 600]);
subplot(2,1,1); loglog(lam, S'); grid on; xlim([15 150]); legend(arrayfun(@(x) sprintf('%g \\mum', x), dz, 'UniformOutput', false), 'Location', 'best');
ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]'); title('Dipole 240 MeV, on axis: tout step sweep');
subplot(2,1,2); loglog(lam, abs(S(1:end-1,:)./S(end,:) - 1)'); grid on; xlim([15 150]); yline(0.01, 'k:');
legend(arrayfun(@(x) sprintf('%g \\mum vs 1.875 \\mum', x), dz(1:end-1), 'UniformOutput', false), 'Location', 'best');
xlabel('\lambda [nm]'); ylabel('|S / S_{1.875\mum} - 1|');
set(findall(fg, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(fg, fullfile(this_dir, 'figures', 'tout_sweep.png'), 'Resolution', 150); savefig(fg, fullfile(this_dir, 'figures', 'tout_sweep.fig'));
