%% make_sec5_figures
% Paper figures for Section 5 (SLAC bunch-length measurement). Nothing is recomputed.
%   fig_slac_scan         single column: energy radiated in 50-300 nm by a single electron on the reference
%                         trajectory along the vertical axis of the observation plane (z = 5.9 m), undulator and
%                         dipole separately, x and y polarization (undulator turned over, dipole bends to +y)
%   fig_slac_calibration  two columns, at (0, -2 mm, 5.9 m) behind an x polarizer, 100 pC:
%                         (a) energy in D1 (50-280 nm) and D2 (300-1000 nm) vs sigma_z, +-1 jackknife sd, incoherent level
%                         (b) spectra for a few sigma_z, with the two detector bands and the incoherent level
%                         (the sensitivity panel was removed on 2026-11-06; data still in fig_slac_calibration.txt)
%   fig_slac_polarizer    single column: spectra at sigma_z = 10 nm without polarizer (undulator alone, dipole alone,
%                         both together; x + y + z) and the total behind the x polarizer
% Sources (read only), matlab/radiation_multi_electron/E7GeV_und_dipole/results/:
%   scan_lines_y.txt          (scan_lines.m; columns 2-5 = turned-over undulator: Wund_x Wund_y Wdip_x Wdip_y [J/sr])
%   large_sigz_summary.mat    (analyze_large_sigz.m; expected value over the longitudinal positions, 5e4 macro-particles)
%   large_sigz_spectra.txt    (incoherent spectrum of the real beam, last column)
%   large_sigz.mat            (run_large_sigz.m; per-chunk sums of the fields shifted to sigma_z = 0, all parts and components)
% Outputs: data/fig_slac_*.txt, data/sec5_slac_figures.mat, figures/fig_slac_*.{pdf,png,fig}
% Run: matlab -batch "make_sec5_figures"

this_dir = fileparts(mfilename('fullpath'));
src = fullfile(this_dir, '..', '..', 'matlab', 'radiation_multi_electron', 'E7GeV_und_dipole', 'results');
data_dir = fullfile(this_dir, 'data');    if ~exist(data_dir, 'dir'); mkdir(data_dir); end
fig_dir  = fullfile(this_dir, 'figures'); if ~exist(fig_dir, 'dir');  mkdir(fig_dir);  end
unused_dir = fullfile(this_dir, 'figures_unused'); if ~exist(unused_dir, 'dir'); mkdir(unused_dir); end  % fig_slac_scan not used in the paper (2026-11-06)
st.fs = 8; st.fs_tick = 7; st.fname = 'Times New Roman';
c1 = [0.000 0.447 0.741]; c2 = [0.850 0.325 0.098]; c_ref = [0.3 0.3 0.3];
y_obs = -2.0;                                           % chosen observation point [mm]

%% data: single-electron scan along y
Y = readmatrix(fullfile(src, 'scan_lines_y.txt'), 'FileType', 'text', 'CommentStyle', '%');
y = Y(:, 1); Wux = Y(:, 2); Wuy = Y(:, 3); Wdx = Y(:, 4); Wdy = Y(:, 5);
i0 = find(abs(y - y_obs) < 1e-9);
fprintf('scan at y = %.1f mm: und x %.3e, und y %.3e, dip x %.3e, dip y %.3e J/sr\n', y(i0), Wux(i0), Wuy(i0), Wdx(i0), Wdy(i0));
fid = fopen(fullfile(data_dir, 'fig_slac_scan.txt'), 'w');
fprintf(fid, '%% Single electron on the reference trajectory, 7 GeV, undulator (turned over) + dipole (bends to +y), observation plane z = 5.9 m, x = 0.\n');
fprintf(fid, '%% Energy radiated in 50-300 nm [J/sr], per sr of the sphere centred at the origin. Source: E7GeV_und_dipole/results/scan_lines_y.txt (sm columns).\n');
fprintf(fid, '%% Columns: 1 y_mm  2 undulator_x  3 undulator_y  4 dipole_x  5 dipole_y\n');
fprintf(fid, '%8.3f %.6e %.6e %.6e %.6e\n', [y Wux Wuy Wdx Wdy]');
fclose(fid);

%% data: calibration curves and spectra
L = load(fullfile(src, 'large_sigz_summary.mat'));
assert(strcmp(L.bn{1}, 'D1 50-280 nm') && strcmp(L.bn{2}, 'D2 300-1000 nm'), 'unexpected band order');
sz = 1e9*L.sz(:); W = L.W(1:2, :)'; Wsd = L.Wsd(1:2, :)'; Winc = L.Winc(1:2)';
sens = 100*L.sens(1:2, :)'; sens_sd = 100*L.sens_sd(1:2, :)';   % stored as fraction per nm -> % per nm
lam = 1e9*L.lam(:);
assert(abs(sens(1, 1) - 17.052) < 0.01, 'sensitivity units: expected 17.05 %%/nm at 10 nm (large_sigz_summary.txt)');
SP = readmatrix(fullfile(src, 'large_sigz_spectra.txt'), 'FileType', 'text', 'CommentStyle', '%');
assert(max(abs(SP(:, 1) - lam)./lam) < 1e-10, 'lambda grids differ');
S_inc = SP(:, end);
sz_show = [10 30 100 300 500];
Ssel = zeros(numel(lam), numel(sz_show));
for k = 1:numel(sz_show); Ssel(:, k) = L.Sall(abs(sz - sz_show(k)) < 1e-6, :)'; end
fid = fopen(fullfile(data_dir, 'fig_slac_calibration.txt'), 'w');
fprintf(fid, '%% (0, -2 mm, 5.9 m), x polarizer, 100 pC, 5e4 macro-particles, expected value over the longitudinal positions (E7GeV_und_dipole/analyze_large_sigz.m).\n');
fprintf(fid, '%% W = energy in the band [J/sr], sd = jackknife standard deviation (25 chunks), sens = -dlnW/dsigma_z [%%/nm]. Incoherent level: D1 %.4e, D2 %.4e J/sr.\n', Winc);
fprintf(fid, '%% Columns: 1 sigma_z_nm  2 W_D1(50-280nm)  3 sd_D1  4 sens_D1  5 sens_sd_D1  6 W_D2(300-1000nm)  7 sd_D2  8 sens_D2  9 sens_sd_D2\n');
fprintf(fid, '%6.1f %.6e %.3e %8.4f %7.4f %.6e %.3e %8.4f %7.4f\n', [sz W(:,1) Wsd(:,1) sens(:,1) sens_sd(:,1) W(:,2) Wsd(:,2) sens(:,2) sens_sd(:,2)]');
fclose(fid);
fid = fopen(fullfile(data_dir, 'fig_slac_spectra.txt'), 'w');
fprintf(fid, '%% (0, -2 mm, 5.9 m), x polarizer, 100 pC, expected value. d2W/domega dOmega [J/(rad/s)/sr].\n');
fprintf(fid, '%% Columns: 1 lambda_nm  then sigma_z = %s nm  last: incoherent level of the real beam\n', num2str(sz_show));
fprintf(fid, [repmat('%.6e ', 1, numel(sz_show) + 1) '%.6e\n'], [lam Ssel S_inc]');
fclose(fid);
save(fullfile(data_dir, 'sec5_slac_figures.mat'), 'y', 'Wux', 'Wuy', 'Wdx', 'Wdy', 'sz', 'W', 'Wsd', 'sens', 'sens_sd', 'Winc', 'lam', 'Ssel', 'sz_show', 'S_inc');

%% Figure: single-electron scan along y
fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 6.0], 'Color', 'w');
tl = tiledlayout(fig, 1, 1, 'Padding', 'compact'); ax = nexttile(tl); hold(ax, 'on');
h1 = semilogy(ax, y, Wux, '-', 'Color', c1, 'LineWidth', 0.9);
h2 = semilogy(ax, y, Wuy, '--', 'Color', c1, 'LineWidth', 0.9);
h3 = semilogy(ax, y, Wdx, '-', 'Color', c2, 'LineWidth', 0.9);
h4 = semilogy(ax, y, Wdy, '--', 'Color', c2, 'LineWidth', 0.9);
xline(ax, y_obs, ':', 'Color', c_ref, 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(ax, 'YScale', 'log', 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', 'Box', 'on', ...
    'XLim', [-10 10], 'XTick', -10:2:10, 'YLim', [1e-26 1e-8], 'YTick', 10.^(-26:4:-10), 'YMinorTick', 'off', 'TickLength', [0.015 0.015]);
xlabel(ax, 'y [mm]', 'FontSize', st.fs);
ylabel(ax, 'dW/d\Omega, 50-300 nm [J sr^{-1}]', 'FontSize', st.fs);
lg = legend(ax, [h1 h2 h3 h4], {'undulator, x', 'undulator, y', 'dipole, x', 'dipole, y'}, 'Location', 'north', ...
    'NumColumns', 2, 'Box', 'off', 'FontName', st.fname, 'FontSize', st.fs_tick);
lg.ItemTokenSize = [14 6];
text(ax, y_obs - 0.3, 3e-25, 'observer', 'HorizontalAlignment', 'right', 'FontName', st.fname, 'FontSize', st.fs_tick - 1, 'Color', c_ref);
text(ax, 9.6, 3e-25, 'bend \rightarrow +y', 'HorizontalAlignment', 'right', 'FontName', st.fname, 'FontSize', st.fs_tick - 1, 'Color', c_ref);
export(fig, fullfile(unused_dir, 'fig_slac_scan'));

%% Figure: calibration, EXPECTED VALUE version (not used since 2026-11-07, written to figures_unused/)
fig = figure('Units', 'centimeters', 'Position', [2 2 17.5 6.5], 'Color', 'w');
tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
cols = [c1; c2]; lab = {'D1, 50-280 nm', 'D2, 300-1000 nm'};
% (a) energy vs sigma_z
ax = nexttile(tl); hold(ax, 'on'); h = gobjects(1, 2);
for b = 1:2
    lo = max(W(:, b) - Wsd(:, b), W(:, b)*1e-3); hi = W(:, b) + Wsd(:, b);
    fill(ax, [sz; flipud(sz)], [lo; flipud(hi)], cols(b, :), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    h(b) = plot(ax, sz, W(:, b), '-', 'Color', cols(b, :), 'LineWidth', 0.9);
    yline(ax, Winc(b), ':', 'Color', cols(b, :), 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
set(ax, 'YScale', 'log', 'XLim', [0 800], 'YLim', [1e-8 1e3], 'YTick', 10.^(-8:2:2), 'YMinorTick', 'off');
xlabel(ax, '\sigma_z [nm]'); ylabel(ax, 'dW/d\Omega in band [J sr^{-1}]');
lg = legend(ax, h, lab, 'Location', 'northeast', 'Box', 'off'); lg.ItemTokenSize = [14 6];
style(ax, st, '(a)');
% (b) spectra (the sensitivity panel was removed on 2026-11-06; three-panel version in figures_unused/)

ax = nexttile(tl); hold(ax, 'on');
yl = [1e-28 1e-13];
for b = 1:2
    bl = 1e9*L.bands(b, :);
    fill(ax, [bl(1) bl(2) bl(2) bl(1)], [yl(1) yl(1) yl(2) yl(2)], cols(b, :), 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
cs = parula(numel(sz_show) + 1); hs = gobjects(1, numel(sz_show) + 1);
for k = 1:numel(sz_show)
    hs(k) = plot(ax, lam, Ssel(:, k), '-', 'Color', cs(k, :), 'LineWidth', 0.8);
end
hs(end) = plot(ax, lam, S_inc, 'k:', 'LineWidth', 1.0);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [50 1000], 'XTick', [50 100 200 300 500 1000], 'YLim', yl, 'YTick', 10.^(-28:2:-14), 'YMinorTick', 'off');
xlabel(ax, '\lambda [nm]'); ylabel(ax, 'd^2W/d\omegad\Omega [J s sr^{-1}]');
lg = legend(ax, hs, [arrayfun(@(s) sprintf('\\sigma_z = %g nm', s), sz_show, 'UniformOutput', false), {'incoherent'}], ...
    'Location', 'south', 'NumColumns', 2, 'Box', 'off'); lg.ItemTokenSize = [12 6];
text(ax, 118, 3e-14, 'D1', 'Color', c1, 'FontName', st.fname, 'FontSize', st.fs_tick, 'HorizontalAlignment', 'center');
text(ax, 550, 3e-14, 'D2', 'Color', c2, 'FontName', st.fname, 'FontSize', st.fs_tick, 'HorizontalAlignment', 'center');
style(ax, st, '(b)', [0.04 0.10]);
export(fig, fullfile(unused_dir, 'fig_slac_calibration_expected_value'));   % replaced by the direct calculation (2026-11-07)

%% Single-column alternative: panel (a) only
fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 6.0], 'Color', 'w');
tl = tiledlayout(fig, 1, 1, 'Padding', 'compact'); ax = nexttile(tl); hold(ax, 'on'); h = gobjects(1, 2);
for b = 1:2
    lo = max(W(:, b) - Wsd(:, b), W(:, b)*1e-3); hi = W(:, b) + Wsd(:, b);
    fill(ax, [sz; flipud(sz)], [lo; flipud(hi)], cols(b, :), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    h(b) = plot(ax, sz, W(:, b), '-', 'Color', cols(b, :), 'LineWidth', 0.9);
    yline(ax, Winc(b), ':', 'Color', cols(b, :), 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
set(ax, 'YScale', 'log', 'XLim', [0 800], 'YLim', [1e-8 1e3], 'YTick', 10.^(-8:2:2), 'YMinorTick', 'off');
xlabel(ax, '\sigma_z [nm]'); ylabel(ax, 'dW/d\Omega in band [J sr^{-1}]');
lg = legend(ax, h, lab, 'Location', 'northeast', 'Box', 'off'); lg.ItemTokenSize = [14 6];
style(ax, st, '');
export(fig, fullfile(unused_dir, 'fig_slac_calibration_single_expected_value'));

%% Figure: calibration, DIRECT calculation (user, 2026-11-07)
% Source: E7GeV_und_dipole/results/direct_sigz/summary.mat (analyze_direct_sigz.m): 5e4 macro-particles, a new GPT run at
% every sigma_z (10 nm = the GPT run of run_large_sigz.m), raw S = pref |sum_j E_j|^2 of the x component (und + dip),
% nothing subtracted; jackknife over 25 chunks. Reference lines: macro-particle noise floor pref sum_j|E_j|^2 (mean over
% sigma_z) and the real incoherent level = noise floor / nmacro.
Dd = load(fullfile(src, 'direct_sigz', 'summary.mat'));
szd = 1e9*Dd.sz(:); Wd = Dd.W'; Wdsd = Dd.Wsd'; Wnd = mean(Dd.Wn, 2)'; lamd = 1e9*Dd.lam(:);
Snd = mean(Dd.Sn, 1)'; Sincd = Snd/Dd.nmac;
fprintf('direct calibration: sigma_z = %s nm; noise floor spread over sigma_z: D1 %.1e, D2 %.1e (rel.)\n', mat2str(szd'), ...
    (max(Dd.Wn(1, :)) - min(Dd.Wn(1, :)))/Wnd(1), (max(Dd.Wn(2, :)) - min(Dd.Wn(2, :)))/Wnd(2));
assert(abs(Wd(1, 1)/82.073 - 1) < 1e-4 && abs(Wd(1, 2)/0.14574 - 1) < 1e-4, 'direct W at 10 nm differs from summary.txt');
szd_show = [10 30 50 100 200];
Sdsel = zeros(numel(lamd), numel(szd_show));
for k = 1:numel(szd_show); Sdsel(:, k) = Dd.S(abs(szd - szd_show(k)) < 1e-6, :)'; end
fid = fopen(fullfile(data_dir, 'fig_slac_calibration_direct.txt'), 'w');
fprintf(fid, '%% (0, -2 mm, 5.9 m), x polarizer, 100 pC, DIRECT calculation: 5e4 macro-particles, new GPT run per sigma_z (E7GeV_und_dipole/results/direct_sigz/summary.mat).\n');
fprintf(fid, '%% W = raw energy in the band [J/sr] (macro-particle self term included, nothing subtracted), sd = jackknife (25 chunks), noise = pref sum_j|E_j|^2 in the band.\n');
fprintf(fid, '%% Mean noise floor: D1 %.4e, D2 %.4e J/sr; real incoherent level = noise / nmacro (%.0f): D1 %.4e, D2 %.4e J/sr.\n', Wnd, Dd.nmac, Wnd/Dd.nmac);
fprintf(fid, '%% Columns: 1 sigma_z_nm  2 W_D1(50-280nm)  3 sd_D1  4 noise_D1  5 W_D2(300-1000nm)  6 sd_D2  7 noise_D2\n');
fprintf(fid, '%6.1f %.6e %.3e %.4e %.6e %.3e %.4e\n', [szd Wd(:,1) Wdsd(:,1) Dd.Wn(1,:)' Wd(:,2) Wdsd(:,2) Dd.Wn(2,:)']');
fclose(fid);
fid = fopen(fullfile(data_dir, 'fig_slac_spectra_direct.txt'), 'w');
fprintf(fid, '%% (0, -2 mm, 5.9 m), x polarizer, 100 pC, DIRECT calculation (raw pref |sum E|^2). d2W/domega dOmega [J/(rad/s)/sr].\n');
fprintf(fid, '%% Columns: 1 lambda_nm  then sigma_z = %s nm  then macro-particle noise floor (mean over sigma_z), real incoherent level\n', num2str(szd_show));
fprintf(fid, [repmat('%.6e ', 1, numel(szd_show) + 2) '%.6e\n'], [lamd Sdsel Snd Sincd]');
fclose(fid);
save(fullfile(data_dir, 'fig_slac_calibration_direct.mat'), 'szd', 'Wd', 'Wdsd', 'Wnd', 'lamd', 'Sdsel', 'szd_show', 'Snd', 'Sincd');

fig = figure('Units', 'centimeters', 'Position', [2 2 17.5 6.5], 'Color', 'w');
tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% (a) energy vs sigma_z
ax = nexttile(tl); hold(ax, 'on'); h = gobjects(1, 4);
for b = 1:2
    h(b) = errorbar(ax, szd, Wd(:, b), Wdsd(:, b), 'o-', 'Color', cols(b, :), 'MarkerFaceColor', cols(b, :), 'MarkerSize', 3, 'LineWidth', 0.9, 'CapSize', 3);
    yline(ax, Wnd(b), '--', 'Color', cols(b, :), 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline(ax, Wnd(b)/Dd.nmac, ':', 'Color', cols(b, :), 'LineWidth', 0.9, 'HandleVisibility', 'off');
end
h(3) = plot(ax, nan, nan, 'k--', 'LineWidth', 0.8); h(4) = plot(ax, nan, nan, 'k:', 'LineWidth', 0.9);
set(ax, 'YScale', 'log', 'XLim', [0 210], 'XTick', 0:50:200, 'YLim', [1e-8 1e4], 'YTick', 10.^(-8:2:4), 'YMinorTick', 'off');
xlabel(ax, '\sigma_z [nm]'); ylabel(ax, 'dW/d\Omega in band [J sr^{-1}]');
lg = legend(ax, h, [lab, {'macro-particle noise floor', 'incoherent (100 pC)'}], 'Location', 'northeast', 'NumColumns', 2, 'Box', 'off'); lg.ItemTokenSize = [12 6];
style(ax, st, '(a)', [0.04 0.10]);
% (b) spectra
ax = nexttile(tl); hold(ax, 'on');
yl = [1e-28 1e-13];
for b = 1:2
    bl = 1e9*L.bands(b, :);
    fill(ax, [bl(1) bl(2) bl(2) bl(1)], [yl(1) yl(1) yl(2) yl(2)], cols(b, :), 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
cs = parula(numel(szd_show) + 1); hs = gobjects(1, numel(szd_show) + 2);
for k = 1:numel(szd_show)
    hs(k) = plot(ax, lamd, Sdsel(:, k), '-', 'Color', cs(k, :), 'LineWidth', 0.8);
end
hs(end-1) = plot(ax, lamd, Snd, 'k--', 'LineWidth', 1.0);
hs(end) = plot(ax, lamd, Sincd, 'k:', 'LineWidth', 1.0);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [50 1000], 'XTick', [50 100 200 300 500 1000], 'YLim', yl, 'YTick', 10.^(-28:2:-14), 'YMinorTick', 'off');
xlabel(ax, '\lambda [nm]'); ylabel(ax, 'd^2W/d\omegad\Omega [J s sr^{-1}]');
lg = legend(ax, hs, [arrayfun(@(s) sprintf('\\sigma_z = %g nm', s), szd_show, 'UniformOutput', false), {'noise floor', 'incoherent'}], ...
    'Location', 'south', 'NumColumns', 2, 'Box', 'off'); lg.ItemTokenSize = [12 6];
text(ax, 118, 3e-14, 'D1', 'Color', c1, 'FontName', st.fname, 'FontSize', st.fs_tick, 'HorizontalAlignment', 'center');
text(ax, 550, 3e-14, 'D2', 'Color', c2, 'FontName', st.fname, 'FontSize', st.fs_tick, 'HorizontalAlignment', 'center');
style(ax, st, '(b)', [0.04 0.10]);
export(fig, fullfile(fig_dir, 'fig_slac_calibration'));

%% Figure: effect of the polarizer at sigma_z = 10 nm
% DIRECT calculation (user, 2026-11-06): one realization of the 5e4 macro-particles of the GPT run at sigma_z = 10 nm,
%   S = pref |sum_j E_j|^2, fields carry q*nmacro (no averaging over the longitudinal positions).
% R.Er(chunk, s, part, comp, omega) = sum_j E_j exp(i w (sz_real(s) - sz_gpt) X3_j / v_zj); for sz_real = sz_gpt = 10 nm
% the phase is exactly 1, i.e. the fields of the GPT run itself.
% undulator: part 1; dipole: part 2; total: parts 1 + 2. Without polarizer: |E_x|^2 + |E_y|^2 + |E_z|^2; x polarizer: |E_x|^2 of the total.
% The result contains the macro-particle self term pref*sum_j|E_j|^2 (nmacro times the real incoherent level); it is
% drawn as the noise level for the x-polarized total. Jackknife over the 25 chunks (fields x 25/24 to keep 100 pC).
% The previous expected-value version is in figures_unused/fig_slac_polarizer_expected_value.*.
c0 = 299792458; eps0 = 8.85418781e-12; sz_p = 10e-9;
R = load(fullfile(src, 'large_sigz.mat')); R = R.R;
dn = find(R.done); nch = numel(dn); assert(nch == numel(R.done), 'not all chunks done');
i10 = find(abs(R.sz_real - sz_p) < 1e-12); assert(numel(i10) == 1 && abs(R.sz_gpt - sz_p) < 1e-12, 'sigma_z = 10 nm realization not found');
lamP = 1e9*R.lambda(:); omP = R.omega(:); prefP = eps0*c0*norm(R.obs)^2/pi;
Er = @(p, c, m) sum(squeeze(R.Er(dn(m), i10, p, c, :)), 1).';       % summed over the chunks m, column over omega
Sdir = @(Ef, f) prefP*(abs(f*Ef(1)).^2 + abs(f*Ef(2)).^2 + abs(f*Ef(3)).^2);
Sfun = @(m, f) [Sdir(@(c) Er(1, c, m), f), Sdir(@(c) Er(2, c, m), f), Sdir(@(c) Er(1, c, m) + Er(2, c, m), f), ...
                prefP*abs(f*(Er(1, 1, m) + Er(2, 1, m))).^2];
P = Sfun(1:nch, 1);                                                  % columns: undulator, dipole, total (no polarizer), total x
J = zeros(nch, numel(lamP), 4);
for c = 1:nch; J(c, :, :) = Sfun(setdiff(1:nch, c), nch/(nch - 1)); end
Psd = squeeze(sqrt((nch - 1)/nch*sum((J - mean(J, 1)).^2, 1)));
Snoise = prefP*sum(squeeze(R.I(dn, 3, 1, :)), 1).';                  % macro-particle self term, x component of und + dip
WbP = @(S, a, b) -trapz(omP(lamP >= a & lamP <= b), S(lamP >= a & lamP <= b));
% checks against the earlier direct sums of the same macro-particles
chk = [WbP(P(:, 4), 50, 280)/8.2073e1 - 1, WbP(P(:, 4), 300, 1000)/1.4574e-1 - 1];   % large_sigz_summary.txt, CHECK 1, realisation at 10 nm
fprintf('polarizer figure (direct): x-polarized total vs CHECK 1 realisation at 10 nm: D1 %+.2e, D2 %+.2e\n', chk);
assert(all(abs(chk) < 1e-4), 'direct x-polarized total does not reproduce the stored realisation');
chk2 = [WbP(P(:, 1), 50, 300)/2.1984e3, WbP(P(:, 2), 50, 300)/4.6304e2, WbP(P(:, 3), 50, 300)/1.5497e3] - 1;   % W_vs_sigz_y-2.0mm.txt (combined.mat, 300 points)
fprintf('  unpolarized, 50-300 nm vs combined.mat (other lambda grid): undulator %+.2e, dipole %+.2e, total %+.2e\n', chk2);
assert(all(abs(chk2) < 5e-3), 'direct unpolarized spectra do not reproduce combined.mat');
bandsP = [50 280; 300 1000];
for b = 1:2
    Wp = arrayfun(@(q) WbP(P(:, q), bandsP(b, 1), bandsP(b, 2)), 1:4);
    fprintf('  band %g-%g nm, W [J/sr]: undulator %.4g, dipole %.4g, total %.4g, total x %.4g (x / total %.3g); x noise level %.3g (%.2g of total x)\n', ...
        bandsP(b, :), Wp, Wp(4)/Wp(3), WbP(Snoise, bandsP(b, 1), bandsP(b, 2)), WbP(Snoise, bandsP(b, 1), bandsP(b, 2))/Wp(4));
end
fprintf('  spectral noise / total x: %.2g at 50 nm, %.2g at 300 nm, %.2g at 1000 nm (max %.2g)\n', Snoise(1)/P(1, 4), ...
    interp1(lamP, Snoise./P(:, 4), 300), Snoise(end)/P(end, 4), max(Snoise./P(:, 4)));
fid = fopen(fullfile(data_dir, 'fig_slac_polarizer.txt'), 'w');
fprintf(fid, '%% (0, -2 mm, 5.9 m), sigma_z = 10 nm, 100 pC, DIRECT sum of one realization: 5e4 macro-particles of the GPT run (E7GeV_und_dipole/results/large_sigz.mat, Er at 10 nm).\n');
fprintf(fid, '%% d2W/domega dOmega [J/(rad/s)/sr]. No polarizer: x + y + z of the undulator alone, the dipole alone, and both together; x polarizer: x component of both together.\n');
fprintf(fid, '%% Contains the macro-particle self term; noise level = pref sum_j|E_j|^2 of the x-polarized total (nmacro = %.0f times the real incoherent level).\n', R.nmacro);
fprintf(fid, '%% sd = jackknife standard deviation (25 chunks). Columns: 1 lambda_nm  2 undulator  3 dipole  4 total  5 total_x  6 sd_undulator  7 sd_dipole  8 sd_total  9 sd_total_x  10 noise_level_x\n');
fprintf(fid, '%.6f %.6e %.6e %.6e %.6e %.3e %.3e %.3e %.3e %.6e\n', [lamP P Psd Snoise]');
fclose(fid);
save(fullfile(data_dir, 'fig_slac_polarizer.mat'), 'lamP', 'P', 'Psd', 'Snoise', 'sz_p');

fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 6.5], 'Color', 'w');
tl = tiledlayout(fig, 1, 1, 'Padding', 'compact'); ax = nexttile(tl); hold(ax, 'on');
yl = [1e-26 1e-12];
for b = 1:2
    fill(ax, [bandsP(b, 1) bandsP(b, 2) bandsP(b, 2) bandsP(b, 1)], [yl(1) yl(1) yl(2) yl(2)], cols(b, :), 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
hp = gobjects(1, 5);
hp(1) = plot(ax, lamP, P(:, 1), '-', 'Color', c1, 'LineWidth', 0.8);
hp(2) = plot(ax, lamP, P(:, 2), '-', 'Color', c2, 'LineWidth', 0.8);
hp(3) = plot(ax, lamP, P(:, 3), 'k-', 'LineWidth', 1.0);
hp(4) = plot(ax, lamP, P(:, 4), '-', 'Color', [0.466 0.674 0.188], 'LineWidth', 1.2);
hp(5) = plot(ax, lamP, Snoise, ':', 'Color', [0.466 0.674 0.188], 'LineWidth', 1.0);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [50 1000], 'XTick', [50 100 200 300 500 1000], 'YLim', yl, 'YTick', 10.^(-26:2:-12), 'YMinorTick', 'off');
xlabel(ax, '\lambda [nm]'); ylabel(ax, 'd^2W/d\omegad\Omega [J s sr^{-1}]');
lg = legend(ax, hp, {'undulator', 'dipole', 'undulator + dipole', 'undulator + dipole, x polarizer', 'macro-particle noise level, x'}, 'Location', 'southwest', 'Box', 'off');
lg.ItemTokenSize = [14 6];
lg.Title.String = '\sigma_z = 10 nm'; lg.Title.FontWeight = 'normal';
text(ax, 118, 4e-13, 'D1', 'Color', c1, 'FontName', st.fname, 'FontSize', st.fs_tick, 'HorizontalAlignment', 'center');
text(ax, 550, 4e-13, 'D2', 'Color', c2, 'FontName', st.fname, 'FontSize', st.fs_tick, 'HorizontalAlignment', 'center');
style(ax, st, '');
lg.Title.FontSize = st.fs_tick - 1; lg.Title.FontName = st.fname;
export(fig, fullfile(fig_dir, 'fig_slac_polarizer'));
fprintf('Figures written to %s\n', fig_dir);

function style(ax, st, letter, pos)
if nargin < 4; pos = [0.04 0.95]; end
set(ax, 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', 'Box', 'on', 'TickLength', [0.015 0.015], 'Layer', 'top');
ax.XLabel.FontSize = st.fs; ax.YLabel.FontSize = st.fs;
lg = ax.Legend; if ~isempty(lg); lg.FontName = st.fname; lg.FontSize = st.fs_tick - 1; end
text(ax, pos(1), pos(2), letter, 'Units', 'normalized', 'FontName', st.fname, 'FontSize', st.fs, 'VerticalAlignment', 'top');
end

function export(fig, base)
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [base '.png'], 'Resolution', 600);
savefig(fig, [base '.fig']);
close(fig);
end
