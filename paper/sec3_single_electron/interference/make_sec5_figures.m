%% make_sec5_figures
% Paper figures for the undulator + dipole interference (single electron, 240 MeV):
%   fig_interference_spectra_vs_D  (Fig. 5) two panels in one column:
%       (a) on axis, D = 0.5 m, 60-150 nm: |E_x +- E_y|^2/2 and (S_und + S_dip)/2
%           (legend as formulas: the +-45 deg projections are selected from the computed fields, not measured);
%           black triangle on the axis at 70 nm = wavelength of panel (b)
%       (b) on axis, lambda = 70 nm: |E_x +- E_y|^2/2 vs D (MAGNETO) with the geometric prediction
%   fig_interference_spectra, fig_interference_vs_D: the two panels as separate figures, now in figures_unused/
%   fig_interference_delay    (Fig. 6) delay of the dipole light vs D: measured from the fringes and predicted, on axis only
%   fig_interference_map      not in the paper (figures_unused/): cos(dphi) in the (lambda, D) plane, on axis
%
% Data (read only): data/spectra_D0p5m_*.txt, data/delay_vs_D_*.txt, data/sec5_interference_maps.mat,
%                   data/fringe_vs_D_70nm_onaxis.txt (collect_sec5_fringe_vs_D.m)
% Run: matlab -batch "make_sec5_figures"

this_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(this_dir, 'data');
fig_dir  = fullfile(this_dir, 'figures'); if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

st.fs = 8; st.fs_tick = 7; st.fname = 'Times New Roman';
c1 = [0.000 0.447 0.741];  c2 = [0.850 0.325 0.098];  c_ref = [0.3 0.3 0.3];
rd = @(f) readmatrix(fullfile(data_dir, f), 'FileType', 'text', 'CommentStyle', '%');
On = rd('spectra_D0p5m_onaxis.txt');  Th = rd('spectra_D0p5m_th0p5.txt');
% columns: 1 lambda 2 omega 3 S0 4 Sx 5 Sy 6 S(+45) 7 S(-45) 8 S2 9 S3 10 und alone 11 dip alone 12 S0-und-dip
lam = On(:, 1);
win = [60 150];

%% Figure 5: spectra at D = 0.5 m, on axis behind a +-45 deg polariser (log y)
fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 6.0], 'Color', 'w');
tl  = tiledlayout(fig, 1, 1, 'Padding', 'compact'); ax = nexttile(tl);
h3 = loglog(ax, lam, (On(:, 10) + On(:, 11)) / 2, '-', 'Color', c_ref, 'LineWidth', 1.2); hold(ax, 'on');
h1 = loglog(ax, lam, On(:, 6), '-', 'Color', c1, 'LineWidth', 0.6);
h2 = loglog(ax, lam, On(:, 7), '-', 'Color', c2, 'LineWidth', 0.6);
lg = legend(ax, [h1 h2 h3], {'|E_x + E_y|^2/2', '|E_x - E_y|^2/2', '(S_{und} + S_{dip})/2'}, ...
    'Location', 'northwest', 'Box', 'off', 'FontName', st.fname, 'FontSize', st.fs_tick);
lg.ItemTokenSize = [14 6];
lg.Title.String = 'on axis, D = 0.5 m'; lg.Title.FontWeight = 'normal'; lg.Title.FontSize = st.fs_tick;
set(ax, 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', 'Box', 'on', ...
    'XLim', win, 'XTick', [60 70 80 100 120 150], 'XMinorTick', 'off', 'TickLength', [0.015 0.015], ...
    'YLim', [1e-34 1e-28], 'YTick', 10.^(-34:2:-28), 'YMinorTick', 'off');
ylabel(ax, 'd^2W/d\omegad\Omega [J s sr^{-1}]', 'FontSize', st.fs);
xlabel(ax, '\lambda [nm]', 'FontSize', st.fs);
export(fig, fullfile(this_dir, 'figures_unused', 'fig_interference_spectra'));   % superseded by fig_interference_spectra_vs_D (a)

%% Energy at a fixed wavelength vs D: lambda = 70 nm, on axis
F = rd('fringe_vs_D_70nm_onaxis.txt');
% columns: 1 D 2 |Ex+Ey|^2/2 3 |Ex-Ey|^2/2 4 (S_und+S_dip)/2 5 prediction +45 6 prediction -45 7 c*dtau
sc = 1e-33; k = 1:3:size(F, 1);                                      % markers every 3 mm
fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 6.0], 'Color', 'w');
tl  = tiledlayout(fig, 1, 1, 'Padding', 'compact'); ax = nexttile(tl);
h3 = plot(ax, F(:, 1), F(:, 4)/sc, '-', 'Color', c_ref, 'LineWidth', 0.8); hold(ax, 'on');
plot(ax, F(:, 1), F(:, 5)/sc, '-', 'Color', c1, 'LineWidth', 0.6);
plot(ax, F(:, 1), F(:, 6)/sc, '-', 'Color', c2, 'LineWidth', 0.6);
h1 = plot(ax, F(k, 1), F(k, 2)/sc, 'o', 'Color', c1, 'MarkerSize', 2.5, 'LineWidth', 0.5);
h2 = plot(ax, F(k, 1), F(k, 3)/sc, 's', 'Color', c2, 'MarkerSize', 2.5, 'LineWidth', 0.5);
h4 = plot(ax, nan, nan, 'k-', 'LineWidth', 0.6);
set(ax, 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', 'Box', 'on', ...
    'XLim', [0.5 0.7], 'YLim', [0 10], 'YTick', 0:2:8, 'TickLength', [0.015 0.015]);
xlabel(ax, 'D [m]', 'FontSize', st.fs);
ylabel(ax, 'd^2W/d\omegad\Omega [10^{-33} J s sr^{-1}]', 'FontSize', st.fs);
lg = legend(ax, [h1 h2 h4 h3], {'|E_x + E_y|^2/2', '|E_x - E_y|^2/2', 'geometric delay', '(S_{und} + S_{dip})/2'}, ...
    'Location', 'north', 'NumColumns', 2, 'Box', 'off', 'FontName', st.fname, 'FontSize', st.fs_tick - 1);
lg.ItemTokenSize = [12 6];
lg.Title.String = 'on axis, \lambda = 70 nm'; lg.Title.FontWeight = 'normal'; lg.Title.FontSize = st.fs_tick;
export(fig, fullfile(this_dir, 'figures_unused', 'fig_interference_vs_D'));       % superseded by fig_interference_spectra_vs_D (b)

%% Paper figure: (a) spectra at D = 0.5 m and (b) energy at 70 nm vs D, stacked in one column
fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 11.0], 'Color', 'w');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
axa = nexttile(tl);
h3 = loglog(axa, lam, (On(:, 10) + On(:, 11)) / 2, '-', 'Color', c_ref, 'LineWidth', 1.2); hold(axa, 'on');
h1 = loglog(axa, lam, On(:, 6), '-', 'Color', c1, 'LineWidth', 0.6);
h2 = loglog(axa, lam, On(:, 7), '-', 'Color', c2, 'LineWidth', 0.6);
plot(axa, 70, 1.5e-34, 'k^', 'MarkerSize', 3.5, 'MarkerFaceColor', 'k');       % wavelength of panel (b)
lg = legend(axa, [h1 h2 h3], {'|E_x + E_y|^2/2', '|E_x - E_y|^2/2', '(S_{und} + S_{dip})/2'}, ...
    'Location', 'northwest', 'Box', 'off', 'FontName', st.fname, 'FontSize', st.fs_tick);
lg.ItemTokenSize = [14 6];
lg.Title.String = 'on axis, D = 0.5 m'; lg.Title.FontWeight = 'normal'; lg.Title.FontSize = st.fs_tick;
set(axa, 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', 'Box', 'on', ...
    'XLim', win, 'XTick', [60 70 80 100 120 150], 'XMinorTick', 'off', 'TickLength', [0.015 0.015], ...
    'YLim', [1e-34 1e-28], 'YTick', 10.^(-34:2:-28), 'YMinorTick', 'off');
ylabel(axa, 'd^2W/d\omegad\Omega [J s sr^{-1}]', 'FontSize', st.fs);
xlabel(axa, '\lambda [nm]', 'FontSize', st.fs);

axb = nexttile(tl);
k = 1:3:size(F, 1);
h3 = plot(axb, F(:, 1), F(:, 4)/sc, '-', 'Color', c_ref, 'LineWidth', 0.8); hold(axb, 'on');
plot(axb, F(:, 1), F(:, 5)/sc, '-', 'Color', c1, 'LineWidth', 0.6);
plot(axb, F(:, 1), F(:, 6)/sc, '-', 'Color', c2, 'LineWidth', 0.6);
h1 = plot(axb, F(k, 1), F(k, 2)/sc, 'o', 'Color', c1, 'MarkerSize', 2.5, 'LineWidth', 0.5);
h2 = plot(axb, F(k, 1), F(k, 3)/sc, 's', 'Color', c2, 'MarkerSize', 2.5, 'LineWidth', 0.5);
h4 = plot(axb, nan, nan, 'k-', 'LineWidth', 0.6);
set(axb, 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', 'Box', 'on', ...
    'XLim', [0.5 0.7], 'YLim', [0 10], 'YTick', 0:2:8, 'TickLength', [0.015 0.015]);
xlabel(axb, 'D [m]', 'FontSize', st.fs);
ylabel(axb, 'd^2W/d\omegad\Omega [10^{-33} J s sr^{-1}]', 'FontSize', st.fs);
lg = legend(axb, [h1 h2 h4 h3], {'|E_x + E_y|^2/2', '|E_x - E_y|^2/2', 'geometric delay', '(S_{und} + S_{dip})/2'}, ...
    'Location', 'north', 'NumColumns', 2, 'Box', 'off', 'FontName', st.fname, 'FontSize', st.fs_tick - 1);
lg.ItemTokenSize = [12 6];
lg.Title.String = 'on axis, \lambda = 70 nm'; lg.Title.FontWeight = 'normal'; lg.Title.FontSize = st.fs_tick;

text(axa, 0.97, 0.95, '(a)', 'Units', 'normalized', 'FontName', st.fname, 'FontSize', st.fs, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(axb, 0.97, 0.95, '(b)', 'Units', 'normalized', 'FontName', st.fname, 'FontSize', st.fs, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
export(fig, fullfile(fig_dir, 'fig_interference_spectra_vs_D'));

%% Figure 6: delay vs D (on axis; the 0.5 mrad observer is not in the paper)
Don = rd('delay_vs_D_onaxis.txt');
% columns: 1 D 2 c*dtau measured [um] 3 predicted 4 difference
fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 6.0], 'Color', 'w');
tl  = tiledlayout(fig, 1, 1, 'Padding', 'compact'); ax = nexttile(tl);
[Ds, i] = sort(Don(:, 1));
hp1 = plot(ax, Ds, Don(i, 3), '-', 'Color', c1, 'LineWidth', 0.7); hold(ax, 'on');
coarse_on = abs(Don(:, 1) * 10 - round(Don(:, 1) * 10)) < 1e-9;          % every 0.1 m
fine_on = ~coarse_on & mod(round(Don(:, 1) * 1000), 10) == 0;            % every 10 mm in 0.5-0.7 m
hm1 = plot(ax, Don(coarse_on | fine_on, 1), Don(coarse_on | fine_on, 2), 'o', 'Color', c1, 'MarkerSize', 3, 'LineWidth', 0.5);
set(ax, 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', 'Box', 'on', ...
    'XLim', [0 1.6], 'TickLength', [0.015 0.015]);
xlabel(ax, 'D [m]', 'FontSize', st.fs);
ylabel(ax, 'c\Delta\tau relative to D = 0.5 m [\mum]', 'FontSize', st.fs);
lg = legend(ax, [hm1 hp1], {'from fringes (MAGNETO)', 'geometry'}, 'Location', 'northwest', 'Box', 'off', ...
    'FontName', st.fname, 'FontSize', st.fs_tick - 1);
lg.ItemTokenSize = [12 6];
lg.Title.String = 'on axis'; lg.Title.FontWeight = 'normal'; lg.Title.FontSize = st.fs_tick - 1;
export(fig, fullfile(fig_dir, 'fig_interference_delay'));

%% Not in the paper: (lambda, D) map of cos(dphi), on axis
m = load(fullfile(data_dir, 'sec5_interference_maps.mat'));  M = m.maps.onaxis;
selD = M.D_m >= 0.5 - 1e-9 & M.D_m <= 0.7 + 1e-9;
selL = M.lambda_nm >= 100 & M.lambda_nm <= 150;
[Dsel, k] = sort(M.D_m(selD));
Z = M.cos_dphi(selD, selL);  Z = Z(k, :);
fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 6.0], 'Color', 'w');
tl  = tiledlayout(fig, 1, 1, 'Padding', 'compact'); ax = nexttile(tl);
imagesc(ax, M.lambda_nm(selL), Dsel, Z); axis(ax, 'xy');
colormap(ax, redblue(256)); clim(ax, [-1 1]);
cb = colorbar(ax); cb.Label.String = 'cos\Delta\phi'; cb.FontName = st.fname; cb.FontSize = st.fs_tick;
cb.Label.FontSize = st.fs;
set(ax, 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, 'TickDir', 'out', 'Box', 'on');
xlabel(ax, '\lambda [nm]', 'FontSize', st.fs);
ylabel(ax, 'D [m]', 'FontSize', st.fs);
export(fig, fullfile(this_dir, 'figures_unused', 'fig_interference_map'));
fprintf('Figures written to %s\n', fig_dir);

function export(fig, base)
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [base '.png'], 'Resolution', 600);
savefig(fig, [base '.fig']);
close(fig);
end

function c = redblue(n)
% blue (-1) - white (0) - red (+1)
h = floor(n / 2);
up = linspace(0, 1, h)';
c = [[up up ones(h, 1)]; [ones(n - h, 1) flipud(linspace(0, 1, n - h)') flipud(linspace(0, 1, n - h)')]];
end