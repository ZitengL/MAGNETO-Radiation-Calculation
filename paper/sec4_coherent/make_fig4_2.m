%% make_fig4_2
% Paper figure 4.2: cold bunch (N = 2000, sigma_z = 10 nm), dipole, on axis at 11 m.
% Single panel, same style and colours as fig4_1 (make_fig4_1.m):
%   MAGNETO |sum_j E_j|^2, analytical Gaussian form factor S1[N + N(N-1)|F|^2],
%   and the limits N S1 (fully incoherent) and N^2 S1 (fully coherent) in light colours.
% Replaces the two-panel draft (old version in figures_unused/fig4_2_dipole_cold_two_panel); the version
% with the trajectory-shift curve is in figures_unused/fig4_2_dipole_cold_with_traj_shift.
% Data (read only): data/dipole_N2000_onaxis.txt (collect_sec4_coherent.m)
%   columns: 1 lambda  3 MAGNETO  4 S1  5 analytic Gaussian
% Run: matlab -batch "make_fig4_2"

this_dir = fileparts(mfilename('fullpath'));
D = readmatrix(fullfile(this_dir, 'data', 'dipole_N2000_onaxis.txt'), 'FileType', 'text', 'CommentStyle', '%');
N = 2000;
lam = D(:, 1);

fs = 8; fs_tick = 7; fname = 'Times New Roman';
c_inc = [0.62 0.76 0.90];      % light blue, as in fig4_1
c_coh = [0.78 0.88 0.62];      % light green, as in fig4_1

fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 6.5], 'Color', 'w');
tl  = tiledlayout(fig, 1, 1, 'Padding', 'compact'); ax = nexttile(tl);
h5 = loglog(ax, lam, N^2 * D(:, 4), '-', 'Color', c_coh, 'LineWidth', 0.6); hold(ax, 'on');
h4 = loglog(ax, lam, N * D(:, 4), '-', 'Color', c_inc, 'LineWidth', 0.6);
h1 = loglog(ax, lam, D(:, 3), 'k-', 'LineWidth', 0.8);
h2 = loglog(ax, lam, D(:, 5), '--', 'Color', [0.85 0.10 0.10], 'LineWidth', 0.6);
set(ax, 'FontName', fname, 'FontSize', fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', 'Box', 'on', ...
    'XLim', [15 150], 'XTick', [15 20 30 50 100 150], 'XMinorTick', 'on', 'TickLength', [0.015 0.015], ...
    'YLim', [1e-36 1e-22], 'YTick', 10.^(-36:2:-22), 'YMinorTick', 'off');
xlabel(ax, '\lambda [nm]', 'FontSize', fs);
ylabel(ax, 'd^2W/d\omegad\Omega [J s sr^{-1}]', 'FontSize', fs);
lg = legend(ax, [h1 h2 h4 h5], {'MAGNETO', 'Analytical form factor', 'NS_1', 'N^2S_1'}, ...
    'Location', 'north', 'NumColumns', 2, 'Box', 'off', 'FontName', fname, 'FontSize', fs_tick);
lg.ItemTokenSize = [14 6];

base = fullfile(this_dir, 'figures', 'fig4_2_dipole_cold');
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [base '.png'], 'Resolution', 600);
savefig(fig, [base '.fig']);
close(fig);
fprintf('Written %s\n', base);
