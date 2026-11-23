%% make_sec3_single_electron
% Collects the data used in the paper section "single-electron radiation"
% (MAGNETO vs SRW, 240 MeV undulator and dipole) into ./data and draws the
% paper figures (single column, 8.5 cm; undulator in linear-y and log-y
% versions, dipole in log-y) into ./figures.
%
% Sources (read only, nothing is recomputed):
%   MAGNETO  : matlab/radiation_und15p_240MeV_single_electron/results/und240_spectra.mat
%              matlab/radiation_dipole240_single_electron/results/dipole240_onaxis.mat
%   SRW      : srw/undulator_240MeV_six_spectra_15_150nm/spectrum_theta_*.txt
%              srw/dipole_240MeV_onaxis_15_150nm/dipole240_onaxis_spectrum.txt
%
% SRW gives only Ex, Ey, so all comparisons use the x+y part of MAGNETO.
% Run: matlab -batch "make_sec3_single_electron"

this_dir = fileparts(mfilename('fullpath'));
root     = fullfile(this_dir, '..', '..');
srw_root = fullfile(root, 'srw');
data_dir = fullfile(this_dir, 'data');    if ~exist(data_dir, 'dir'); mkdir(data_dir); end
fig_dir  = fullfile(this_dir, 'figures'); if ~exist(fig_dir, 'dir');  mkdir(fig_dir);  end

%% 1. Undulator: MAGNETO and SRW, six angles
a  = load(fullfile(root, 'matlab', 'radiation_und15p_240MeV_single_electron', 'results', 'und240_spectra.mat'));
sp = a.sp;
lam   = sp.lambda_nm(:);
omega = sp.omega(:);
theta_mrad = sp.theta_x_rad * 1e3;
Na = numel(theta_mrad);

nw  = round(0.01 / mean(diff(log(lam))));            % +-1 % in lambda, same as run_und240_spectra.m
env = @(S) movmax(S, 2*nw + 1);

und = struct();
und.lambda_nm = lam; und.omega = omega; und.theta_x_mrad = theta_mrad;
und.MAGNETO_xy    = zeros(numel(lam), Na);
und.MAGNETO_total = zeros(numel(lam), Na);
und.SRW_xy        = zeros(numel(lam), Na);
und.rel_env       = zeros(numel(lam), Na);
und.rel_point     = zeros(numel(lam), Na);
und_stats = zeros(Na, 7);
for i = 1:Na
    xyz  = squeeze(sp.d2W_dwdOmega_xyz(i, :, :));     % 3 x Nw
    M_xy = (xyz(1, :) + xyz(2, :))';
    fn = fullfile(srw_root, 'undulator_240MeV_six_spectra_15_150nm', sprintf('spectrum_theta_%.1fmrad.txt', theta_mrad(i)));
    A  = readmatrix(fn, 'FileType', 'text', 'CommentStyle', '#');
    % consistency: same wavelength nodes, and SRW file's MATLAB column = this MAGNETO result
    assert(max(abs(A(:, 1) - lam) ./ lam) < 1e-14, 'lambda grid mismatch, %s', fn);
    chk = max(abs(A(:, 6) - M_xy)) / max(M_xy);
    assert(chk < 1e-12, 'MAGNETO column in %s differs from the .mat (%.2e)', fn, chk);
    S_xy = A(:, 3);
    E = env(M_xy);
    und.MAGNETO_xy(:, i)    = M_xy;
    und.MAGNETO_total(:, i) = sp.d2W_dwdOmega(i, :)';
    und.SRW_xy(:, i)        = S_xy;
    und.rel_env(:, i)       = (S_xy - M_xy) ./ E;
    und.rel_point(:, i)     = (S_xy - M_xy) ./ M_xy;
    [pk, ip]  = max(M_xy);
    [~, ipS]  = max(S_xy);
    und_stats(i, :) = [theta_mrad(i), lam(ip), lam(ipS), pk, max(abs(S_xy - M_xy)) / pk, ...
                       max(abs(und.rel_env(:, i))), max(abs(und.rel_point(:, i)))];
end
und.harmonics_nm = sp.harmonics_nm;

%% 2. Dipole: MAGNETO and SRW, on axis
b = load(fullfile(root, 'matlab', 'radiation_dipole240_single_electron', 'results', 'dipole240_onaxis.mat'));
d = b.d;
assert(max(abs(d.lambda_nm(:) - lam) ./ lam) < 1e-14, 'dipole and undulator lambda grids differ');
B = readmatrix(fullfile(srw_root, 'dipole_240MeV_onaxis_15_150nm', 'dipole240_onaxis_spectrum.txt'), ...
    'FileType', 'text', 'CommentStyle', '#');
assert(max(abs(B(:, 1) - lam) ./ lam) < 1e-14, 'dipole lambda grid mismatch');
Md_xy = (d.d2W_dwdOmega_xyz(1, :) + d.d2W_dwdOmega_xyz(2, :))';
chk = max(abs(B(:, 6) - Md_xy)) / max(Md_xy);
assert(chk < 1e-12, 'dipole MAGNETO column differs from the .mat (%.2e)', chk);

dip = struct();
dip.lambda_nm     = lam; dip.omega = omega;
dip.MAGNETO_xy    = Md_xy;
dip.MAGNETO_total = d.d2W_dwdOmega(:);
dip.SRW_xy        = B(:, 3);
dip.LW_rad_plus_vel_xy = B(:, 11);    % MATLAB reference with radiation + velocity term (from the SRW benchmark)
dip.rel_point     = (dip.SRW_xy - Md_xy) ./ Md_xy;
dip.rel_point_LW  = (dip.SRW_xy - dip.LW_rad_plus_vel_xy) ./ dip.LW_rad_plus_vel_xy;
[pkd, ipd] = max(Md_xy);

%% 3. Save data (.mat and .txt)
save(fullfile(data_dir, 'sec3_single_electron.mat'), 'und', 'dip', 'und_stats');

for i = 1:Na
    fn = fullfile(data_dir, sprintf('undulator_theta%.1fmrad.txt', theta_mrad(i)));
    fid = fopen(fn, 'w');
    fprintf(fid, '%% 240 MeV single electron, 15-period undulator, observer 11 m*(sin(theta_x),0,cos(theta_x)), theta_x = %.1f mrad\n', theta_mrad(i));
    fprintf(fid, '%% Energy in J/(rad/s)/sr. x+y = sum of the x and y polarisations (SRW has no Ez). Delta = SRW_xy - MAGNETO_xy.\n');
    fprintf(fid, '%% Local envelope = max of MAGNETO_xy within +-1%% in lambda.\n');
    fprintf(fid, '%% Columns: lambda_nm  omega_rad_s  MAGNETO_xy  MAGNETO_total_xyz  SRW_xy  Delta  Delta/local_envelope  Delta/MAGNETO_xy\n');
    fprintf(fid, '%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n', ...
        [lam, omega, und.MAGNETO_xy(:, i), und.MAGNETO_total(:, i), und.SRW_xy(:, i), ...
         und.SRW_xy(:, i) - und.MAGNETO_xy(:, i), und.rel_env(:, i), und.rel_point(:, i)]');
    fclose(fid);
end

fid = fopen(fullfile(data_dir, 'dipole_onaxis.txt'), 'w');
fprintf(fid, '%% 240 MeV single electron, dipole 0.65 T / 0.1 m (Enge fringe, gap 1 cm), observer on axis (0,0,11) m\n');
fprintf(fid, '%% Energy in J/(rad/s)/sr. x+y = sum of the x and y polarisations (SRW has no Ez).\n');
fprintf(fid, '%% LW_rad_plus_vel_xy: MATLAB reference with radiation AND velocity terms of the LW field (computed in the SRW benchmark folder).\n');
fprintf(fid, '%% Columns: lambda_nm  omega_rad_s  MAGNETO_xy  MAGNETO_total_xyz  SRW_xy  LW_rad_plus_vel_xy  (SRW-MAGNETO)/MAGNETO  (SRW-LW_rad_plus_vel)/LW_rad_plus_vel\n');
fprintf(fid, '%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n', ...
    [lam, omega, dip.MAGNETO_xy, dip.MAGNETO_total, dip.SRW_xy, dip.LW_rad_plus_vel_xy, dip.rel_point, dip.rel_point_LW]');
fclose(fid);

fid = fopen(fullfile(data_dir, 'summary_numbers.txt'), 'w');
fprintf(fid, '%% Undulator, MAGNETO vs SRW (x+y), 15-150 nm, 5000 log-spaced points\n');
fprintf(fid, '%% theta_x_mrad  MAGNETO_peak_nm  SRW_peak_nm  MAGNETO_peak_J/(rad/s)/sr  max|Delta|/peak  max|Delta|/local_envelope  max|Delta|/MAGNETO_pointwise\n');
fprintf(fid, '%.1f %.5f %.5f %.6e %.3e %.3e %.3e\n', und_stats');
fprintf(fid, '%%\n%% Dipole on axis, MAGNETO vs SRW (x+y)\n');
fprintf(fid, 'peak at %.2f nm, %.6e J/(rad/s)/sr; S(15 nm)/peak = %.3e\n', lam(ipd), pkd, Md_xy(1) / pkd);
fprintf(fid, 'max |SRW-MAGNETO|/MAGNETO = %.3e (at %.3f nm)\n', max(abs(dip.rel_point)), lam(find(abs(dip.rel_point) == max(abs(dip.rel_point)), 1)));
fprintf(fid, 'max |SRW-LW(rad+vel)|/LW(rad+vel) = %.3e\n', max(abs(dip.rel_point_LW)));
for band = [15 20; 20 30; 30 50; 50 100; 100 150]'
    k = lam >= band(1) & lam <= band(2);
    fprintf(fid, '  %3d-%3d nm: max |SRW-MAGNETO|/MAGNETO = %.3e, max |SRW-LW(rad+vel)|/LW = %.3e\n', ...
        band(1), band(2), max(abs(dip.rel_point(k))), max(abs(dip.rel_point_LW(k))));
end
fclose(fid);
type(fullfile(data_dir, 'summary_numbers.txt'));

%% 4. Paper figures: single column (8.5 cm), spectrum on top (2/3), difference below (1/3)
% Undulator in two versions (linear and log y) for the user to choose; dipole in one.
st.fs = 8; st.fs_tick = 7; st.fname = 'Times New Roman';
c_on  = [0.000 0.447 0.741];     % on axis
c_off = [0.850 0.325 0.098];     % 0.5 mrad
i_on  = find(abs(theta_mrad - 0.0) < 1e-9);
i_off = find(abs(theta_mrad - 0.5) < 1e-9);

% --- undulator ---
for ysc = {'liny', 'logy'}
    [fig, ax1, ax2] = two_panel_figure(st);
    if strcmp(ysc{1}, 'liny')
        sc = 1e-29;
        h1 = semilogx(ax1, lam, und.MAGNETO_xy(:, i_on)  / sc, '-', 'Color', c_on,  'LineWidth', 0.7); hold(ax1, 'on');
        h2 = semilogx(ax1, lam, und.MAGNETO_xy(:, i_off) / sc, '-', 'Color', c_off, 'LineWidth', 0.7);
        h3 = semilogx(ax1, lam, und.SRW_xy(:, i_on)  / sc, 'k--', 'LineWidth', 0.5);
        semilogx(ax1, lam, und.SRW_xy(:, i_off) / sc, 'k--', 'LineWidth', 0.5);
        set(ax1, 'YLim', [0 2], 'YTick', 0:0.5:2);
        ylabel(ax1, 'd^2W/d\omegad\Omega [10^{-29} J s sr^{-1}]');
    else
        h1 = loglog(ax1, lam, und.MAGNETO_xy(:, i_on),  '-', 'Color', c_on,  'LineWidth', 0.6); hold(ax1, 'on');
        h2 = loglog(ax1, lam, und.MAGNETO_xy(:, i_off), '-', 'Color', c_off, 'LineWidth', 0.6);
        h3 = loglog(ax1, lam, und.SRW_xy(:, i_on),  'k--', 'LineWidth', 0.4);
        loglog(ax1, lam, und.SRW_xy(:, i_off), 'k--', 'LineWidth', 0.4);
        set(ax1, 'YLim', [1e-34 1e-27], 'YTick', 10.^(-34:1:-27), 'YMinorTick', 'off');
        ylabel(ax1, 'd^2W/d\omegad\Omega [J s sr^{-1}]');
    end
    lg = legend(ax1, [h1 h2 h3], {'MAGNETO, \theta_x = 0', 'MAGNETO, \theta_x = 0.5 mrad', 'SRW'}, ...
        'Location', 'northwest', 'Box', 'off', 'FontName', st.fname, 'FontSize', st.fs_tick);
    lg.ItemTokenSize = [14 6];
    semilogx(ax2, lam, 100 * und.rel_env(:, i_on),  '-', 'Color', c_on,  'LineWidth', 0.5); hold(ax2, 'on');
    semilogx(ax2, lam, 100 * und.rel_env(:, i_off), '-', 'Color', c_off, 'LineWidth', 0.5);
    set(ax2, 'YLim', [-0.1 0.1], 'YTick', -0.1:0.05:0.1);
    ylabel(ax2, '\DeltaS / envelope [%]');
    finish_figure(fig, ax1, ax2, st, fullfile(fig_dir, ['fig_undulator_MAGNETO_vs_SRW_' ysc{1}]));
end

% --- dipole (on axis), log y ---
[fig, ax1, ax2] = two_panel_figure(st);
h1 = loglog(ax1, lam, dip.MAGNETO_xy, '-', 'Color', c_on, 'LineWidth', 0.7); hold(ax1, 'on');
h2 = loglog(ax1, lam, dip.SRW_xy, 'k--', 'LineWidth', 0.5);
set(ax1, 'YLim', [1e-38 1e-32], 'YTick', 10.^(-38:1:-32), 'YMinorTick', 'off');
ylabel(ax1, 'd^2W/d\omegad\Omega [J s sr^{-1}]');
lg = legend(ax1, [h1 h2], {'MAGNETO', 'SRW'}, 'Location', 'southeast', 'Box', 'off', ...
    'FontName', st.fname, 'FontSize', st.fs_tick);
lg.ItemTokenSize = [14 6];
semilogx(ax2, lam, 100 * dip.rel_point, '-', 'Color', c_on, 'LineWidth', 0.6);
set(ax2, 'YLim', [-0.1 0.4], 'YTick', -0.1:0.1:0.4);
ylabel(ax2, '\DeltaS / S [%]');
finish_figure(fig, ax1, ax2, st, fullfile(fig_dir, 'fig_dipole_MAGNETO_vs_SRW'));
fprintf('Figures written to %s\n', fig_dir);

%% local functions
function [fig, ax1, ax2] = two_panel_figure(st)
% 8.5 cm wide; spectrum spans two rows, difference one row
fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 7.0], 'Color', 'w');
tl  = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(tl, 1, [2 1]);
ax2 = nexttile(tl, 3);
set(ax2, 'XScale', 'log');
set([ax1 ax2], 'FontName', st.fname);
end

function finish_figure(fig, ax1, ax2, st, base)
labels = {'(a)', '(b)'};
axs = [ax1 ax2];
for k = 1:2
    ax = axs(k);
    set(ax, 'FontName', st.fname, 'FontSize', st.fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', ...
        'XLim', [15 150], 'XTick', [15 20 30 50 100 150], 'XTickLabelRotation', 0, ...
        'Box', 'on', 'XMinorTick', 'on', 'TickLength', [0.015 0.015]);
    ax.XLabel.FontSize = st.fs; ax.YLabel.FontSize = st.fs;
    text(ax, 0.97, 0.94, labels{k}, 'Units', 'normalized', 'FontName', st.fname, 'FontSize', st.fs, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
end
ax1.XTickLabel = {};
xlabel(ax2, '\lambda [nm]');
linkaxes([ax1 ax2], 'x');
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [base '.png'], 'Resolution', 600);
savefig(fig, [base '.fig']);
end
