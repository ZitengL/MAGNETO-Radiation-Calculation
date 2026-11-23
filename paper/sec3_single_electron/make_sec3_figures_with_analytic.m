%% make_sec3_figures_with_analytic
% Paper figures for Section 3: spectra only (no difference panels; the
% agreement is quoted as numbers in the text, see README):
%   fig_undulator_theta0_logy           MAGNETO, SRW, analytic, on axis (the only undulator figure in the paper)
%   fig_dipole                          two columns: (a) trajectory and the two lines of sight, (b) on axis and
%                                       (c) arc-midpoint tangent: MAGNETO and SRW (fringe-field dipole), analytic (uniform arc)
% Not in the paper (2026-10-27), written to figures_unused/: fig_undulator_theta0_liny, fig_undulator_theta0p5_<liny|logy>.
%
% Data (read only, nothing is recomputed):
%   data/sec3_single_electron.mat           MAGNETO and SRW (x+y), from make_sec3_single_electron.m
%   analytic/data/sec3_analytic.mat         analytic comparison, from analytic/collect_sec3_analytic.m
%   analytic/data/dipole_onaxis_vs_analytic.txt
%   data/dipole_trajectory.mat              from collect_dipole_trajectory.m
%   matlab/dipole_analytic_comparison/dipole240_tangent/results/compare_dipole_tangent.mat (R_src, obs)
%   srw/dipole_240MeV_onaxis_15_150nm/dipole240_tangent/SRW_tangent_enge_11m.txt
%
% All spectra are shown per sr of the 11 m sphere centred at the origin, as in
% the rest of Section 3. The analytic formulas are defined per sr seen from the
% source (undulator centre / tangent point), so they are multiplied by (|obs|/R_src)^2.
% Run: matlab -batch "make_sec3_figures_with_analytic"

this_dir = fileparts(mfilename('fullpath'));
root     = fullfile(this_dir, '..', '..');
fig_dir  = fullfile(this_dir, 'figures');

a  = load(fullfile(this_dir, 'data', 'sec3_single_electron.mat'));          und = a.und; dip = a.dip;
b  = load(fullfile(this_dir, 'analytic', 'data', 'sec3_analytic.mat'));     an  = b.und;
ct = load(fullfile(root, 'matlab', 'dipole_analytic_comparison', 'dipole240_tangent', 'results', 'compare_dipole_tangent.mat'));
ct = ct.out;
lam = und.lambda_nm(:);
assert(max(abs(b.lam(:) - lam) ./ lam) < 1e-14, 'lambda grids differ');

st.fs = 8; st.fs_tick = 7; st.fname = 'Times New Roman';
c_on  = [0.000 0.447 0.741];
c_off = [0.850 0.325 0.098];
c_an  = [0.55 0.55 0.55];
H = 6.0;                                                    % figure height [cm]

%% Undulator, one figure per angle
cases = {'theta0', 0.0, 'c_11m_theta0', c_on; 'theta0p5', 0.5, 'c_11m_theta0p5', c_off};
for k = 1:size(cases, 1)
    i_th = find(abs(und.theta_x_mrad - cases{k, 2}) < 1e-9);
    c    = an.(cases{k, 3});
    conv = (norm(c.obs) / c.R_src)^2;                      % per sr from the undulator centre -> per sr of the 11 m sphere
    S_an = c.M(:, 3) * conv;
    chk  = max(abs(c.M(:, 6) * conv - und.MAGNETO_total(:, i_th))) / max(und.MAGNETO_total(:, i_th));
    assert(chk < 1e-6, 'normalisation check failed for %s (%.2e)', cases{k, 1}, chk);

    for ysc = {'liny', 'logy'}
        [fig, ax] = one_panel_figure(st, H);
        if strcmp(ysc{1}, 'liny')
            sc = 1e-29;
            h3 = semilogx(ax, lam, S_an / sc, '-', 'Color', c_an, 'LineWidth', 0.9); hold(ax, 'on');
            h1 = semilogx(ax, lam, und.MAGNETO_xy(:, i_th) / sc, '-', 'Color', cases{k, 4}, 'LineWidth', 0.6);
            h2 = semilogx(ax, lam, und.SRW_xy(:, i_th) / sc, 'k--', 'LineWidth', 0.5);
            set(ax, 'YLim', [0 2], 'YTick', 0:0.5:2);
            ylabel(ax, 'd^2W/d\omegad\Omega [10^{-29} J s sr^{-1}]');
            lg = legend(ax, [h1 h2 h3], {'MAGNETO', 'SRW', 'analytic'}, 'Location', 'northwest', ...
                'Box', 'off', 'FontName', st.fname, 'FontSize', st.fs_tick);
        else
            h3 = loglog(ax, lam, S_an, '-', 'Color', c_an, 'LineWidth', 0.7); hold(ax, 'on');
            h1 = loglog(ax, lam, und.MAGNETO_xy(:, i_th), '-', 'Color', cases{k, 4}, 'LineWidth', 0.5);
            h2 = loglog(ax, lam, und.SRW_xy(:, i_th), 'k--', 'LineWidth', 0.4);
            set(ax, 'YLim', [1e-34 1e-26], 'YTick', 10.^(-34:2:-26), 'YMinorTick', 'off');
            ylabel(ax, 'd^2W/d\omegad\Omega [J s sr^{-1}]');
            lg = legend(ax, [h1 h2 h3], {'MAGNETO', 'SRW', 'analytic'}, 'Location', 'north', 'NumColumns', 3, ...
                'Box', 'off', 'FontName', st.fname, 'FontSize', st.fs_tick);   % one row in the empty top decades
        end
        lg.ItemTokenSize = [14 6];
        lg.Title.String = sprintf('\\theta_x = %g mrad', cases{k, 2});
        lg.Title.FontWeight = 'normal'; lg.Title.FontSize = st.fs_tick; lg.Title.FontName = st.fname;
        out_dir = fig_dir; if ~(k == 1 && strcmp(ysc{1}, 'logy')); out_dir = fullfile(this_dir, 'figures_unused'); end
        finish_figure(fig, ax, st, fullfile(out_dir, sprintf('fig_undulator_%s_%s', cases{k, 1}, ysc{1})));
    end
end

%% Dipole on axis: MAGNETO, SRW, analytic (formula assumes a uniform arc; not applicable on axis)
A = readmatrix(fullfile(this_dir, 'analytic', 'data', 'dipole_onaxis_vs_analytic.txt'), 'FileType', 'text', 'CommentStyle', '%');
assert(max(abs(A(:, 1) - lam) ./ lam) < 1e-8, 'dipole on-axis lambda grid differs');
r = dip.MAGNETO_total(:) ./ A(:, 3);                      % per sr from the tangent point -> per sr of the 11 m sphere
conv_d = median(r);
assert(max(abs(r / conv_d - 1)) < 1e-5, 'dipole on-axis normalisation factor not constant');
S_an_d = A(:, 4) * conv_d;

%% Dipole, arc-midpoint tangent: MAGNETO and SRW (fringe-field dipole) vs analytic
o = ct.o11m;
assert(max(abs(ct.lambda_nm(:) - lam) ./ lam) < 1e-14, 'tangent lambda grid differs');
conv_t = (norm(o.obs) / o.R_src)^2;
% SRW (native x+y), same observer, per sr seen from the tangent point like the MAGNETO tangent result
T = readmatrix(fullfile(root, 'srw', 'dipole_240MeV_onaxis_15_150nm', 'dipole240_tangent', 'SRW_tangent_enge_11m.txt'), ...
    'FileType', 'text', 'CommentStyle', '#');
assert(max(abs(T(:, 1) - lam) ./ lam) < 1e-14, 'SRW tangent lambda grid differs');
chk = max(abs(T(:, 6) - o.S_MAGNETO_enge(:)) ./ o.S_MAGNETO_enge(:));
assert(chk < 1e-10, 'MAGNETO column of the SRW tangent file differs from compare_dipole_tangent.mat (%.2e)', chk);
S_srw_t = T(:, 3);
fprintf('dipole tangent: max |SRW_xy/MAGNETO_total - 1| = %.3e, max |SRW_xy/analytic - 1| = %.3e, max |SRW_xy/MATLAB_xy - 1| = %.3e, max |MAGNETO_total/analytic - 1| = %.3e\n', ...
    max(abs(T(:, 19))), max(abs(T(:, 17))), max(abs(T(:, 18))), max(abs(o.ratio_enge_analytic_m1)));

%% Dipole figure (two columns): (a) geometry, (b) on axis, (c) arc-midpoint tangent
P = load(fullfile(this_dir, 'data', 'dipole_trajectory.mat'));
c_tan = c_off;
lam_c = 49.6;                                                % critical wavelength [nm] of the 240 MeV dipole (user, 2026-11-07)
c_lc  = [0.30 0.30 0.30];                                    % dashed vertical line in (b) and (c)
fig = figure('Units', 'centimeters', 'Position', [2 2 17.5 5.6], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% (a) trajectory (y in mm, z in m: not to scale) and the two lines of sight
ax = nexttile(tl); hold(ax, 'on');
zl = [0.85 1.15]; yl = [-2 10];
fill(ax, [0.95 1.05 1.05 0.95], [yl(1) yl(1) yl(2) yl(2)], [0.92 0.92 0.92], 'EdgeColor', 'none');
text(ax, 1.0, 9.2, 'dipole', 'HorizontalAlignment', 'center', 'FontName', st.fname, 'FontSize', st.fs_tick);
plot(ax, zl, [0 0], '--', 'Color', c_on, 'LineWidth', 0.6);                                   % line of sight, on axis
k_t = (P.obs(2, 2) - P.tp(2, 2)) / (P.obs(2, 3) - P.tp(2, 3));                              % towards (0, 0.406788, 11) m
plot(ax, zl, 1e3*(P.tp(2, 2) + k_t*(zl - P.tp(2, 3))), '--', 'Color', c_tan, 'LineWidth', 0.6);
plot(ax, [zl(1); P.z], 1e3*[0; P.y], 'k-', 'LineWidth', 0.8);                              % straight upstream of the output start
cc = [c_on; c_tan];
for i = 1:2                                                  % parts within 1/gamma of each line of sight, where |B| > 1% of B0
    s = P.z >= P.zr(i, 1) & P.z <= P.zr(i, 2) & P.B > 0.01*max(P.B);
    plot(ax, P.z(s), 1e3*P.y(s), '-', 'Color', cc(i, :), 'LineWidth', 2.2);
end
text(ax, 0.86, -1.0, 'to on-axis detector, (b)', 'Color', c_on, 'FontName', st.fname, 'FontSize', st.fs_tick - 1);   % panel letter = spectrum seen there
text(ax, 1.105, 2.2, {'to arc-midpoint', 'detector, (c)'}, 'Color', c_tan, 'FontName', st.fname, 'FontSize', st.fs_tick - 1, 'HorizontalAlignment', 'center');
text(ax, 1.12, 8.3, 'e^-', 'FontName', st.fname, 'FontSize', st.fs_tick);
set(ax, 'XLim', zl, 'YLim', yl, 'XTick', 0.85:0.05:1.15, 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, ...
    'TickDir', 'in', 'Box', 'on', 'TickLength', [0.015 0.015], 'Layer', 'top');
xlabel(ax, 'z [m]', 'FontSize', st.fs); ylabel(ax, 'y [mm]', 'FontSize', st.fs);
text(ax, 0.03, 0.93, '(a)', 'Units', 'normalized', 'FontName', st.fname, 'FontSize', st.fs);

% (b) on axis
ax = nexttile(tl);
h3 = loglog(ax, lam, S_an_d, '-', 'Color', c_an, 'LineWidth', 0.9); hold(ax, 'on');
h1 = loglog(ax, lam, dip.MAGNETO_xy, '-', 'Color', c_on, 'LineWidth', 0.7);
h2 = loglog(ax, lam, dip.SRW_xy, 'k--', 'LineWidth', 0.5);
set(ax, 'YLim', [1e-38 1e-30], 'YTick', 10.^(-38:2:-30), 'YMinorTick', 'off');
ylabel(ax, 'd^2W/d\omegad\Omega [J s sr^{-1}]');
lg = legend(ax, [h1 h2 h3], {'MAGNETO', 'SRW', 'analytic (uniform arc)'}, 'Location', 'southeast', ...
    'Box', 'off', 'FontName', st.fname, 'FontSize', st.fs_tick);
lg.ItemTokenSize = [14 6];
lg.Title.String = 'on axis'; lg.Title.FontWeight = 'normal'; lg.Title.FontSize = st.fs_tick; lg.Title.FontName = st.fname;
style_spectrum(ax, st);
mark_lambda_c(ax, st, lam_c, c_lc);
text(ax, 0.03, 0.93, '(b)', 'Units', 'normalized', 'FontName', st.fname, 'FontSize', st.fs);

% (c) arc-midpoint tangent
ax = nexttile(tl);
h3 = loglog(ax, lam, o.S_analytic * conv_t, '-', 'Color', c_an, 'LineWidth', 1.6); hold(ax, 'on');
h1 = loglog(ax, lam, o.S_MAGNETO_enge * conv_t, '-', 'Color', c_tan, 'LineWidth', 0.7);
h2 = loglog(ax, lam, S_srw_t * conv_t, 'k--', 'LineWidth', 0.5);
set(ax, 'YLim', [1e-33 1e-31], 'YMinorTick', 'on');
lg = legend(ax, [h1 h2 h3], {'MAGNETO', 'SRW', 'analytic (uniform arc)'}, 'Location', 'southeast', ...
    'Box', 'off', 'FontName', st.fname, 'FontSize', st.fs_tick);
lg.ItemTokenSize = [14 6];
lg.Title.String = 'arc-midpoint tangent'; lg.Title.FontWeight = 'normal'; lg.Title.FontSize = st.fs_tick; lg.Title.FontName = st.fname;
style_spectrum(ax, st);
mark_lambda_c(ax, st, lam_c, c_lc);
text(ax, 0.03, 0.93, '(c)', 'Units', 'normalized', 'FontName', st.fname, 'FontSize', st.fs);

base = fullfile(fig_dir, 'fig_dipole');
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [base '.png'], 'Resolution', 600);
savefig(fig, [base '.fig']); close(fig);
fprintf('Figures written to %s\n', fig_dir);

%% local functions
function [fig, ax] = one_panel_figure(st, height_cm)
fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 height_cm], 'Color', 'w');
tl  = tiledlayout(fig, 1, 1, 'Padding', 'compact');
ax  = nexttile(tl);
set(ax, 'FontName', st.fname);
end

function style_spectrum(ax, st)
set(ax, 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', ...
    'XLim', [15 150], 'XTick', [15 20 30 50 100 150], 'XTickLabelRotation', 0, ...
    'Box', 'on', 'XMinorTick', 'on', 'TickLength', [0.015 0.015]);
xlabel(ax, '\lambda [nm]');
ax.XLabel.FontSize = st.fs; ax.YLabel.FontSize = st.fs;
end

function mark_lambda_c(ax, st, lam_c, col)
xl = xline(ax, lam_c, '--', sprintf('\\lambda_c = %.1f nm', lam_c), 'Color', col, 'LineWidth', 0.6, ...
    'HandleVisibility', 'off', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'top', ...
    'LabelHorizontalAlignment', 'left', 'FontName', st.fname, 'FontSize', st.fs_tick - 1);
uistack(xl, 'bottom');
lg = ax.Legend; lg.Box = 'on'; lg.EdgeColor = 'none'; lg.Color = 'w';   % hide the line behind the legend text
end

function finish_figure(fig, ax, st, base)
style_spectrum(ax, st);
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [base '.png'], 'Resolution', 600);
savefig(fig, [base '.fig']);
close(fig);
end
