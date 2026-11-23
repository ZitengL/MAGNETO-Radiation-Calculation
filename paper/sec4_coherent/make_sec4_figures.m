%% make_sec4_figures
% Draft figures for Section 4 (coherent radiation), from ./data only (run collect_sec4_coherent.m and
% collect_sec4_spread_emittance.m first). Draft style: Times New Roman, not yet the final paper size.
%   fig4_1_undulator_cold      draft two-panel version, now written to figures_unused/fig4_1_undulator_cold_two_panel;
%                              the paper version (single panel) is made by make_fig4_1.m
%   fig4_2_dipole_cold         draft two-panel version (spectrum; MAGNETO/analytic - 1 with the +-1 sigma shot-noise band),
%                              now written to figures_unused/fig4_2_dipole_cold_two_panel; paper version: make_fig4_2.m
%   fig4_3_lines_inc_<obs>     incoherent line shapes, normalised to the peak (top: energy spread, bottom: emittance)
%   fig4_4_fwhm_<obs>          FWHM vs sigma_delta and en: incoherent (filled, only where defined), coherent
%                              (open squares, only where r >= 10), simple estimate (dashed)
%   fig4_5_coherence_<obs>     r = int S_coh / int S_inc vs sigma_delta and en, with the form-factor-model prediction
% Only the on-axis versions go to figures/; the 0.5 mrad versions are written to figures_unused/ (not in the paper).
% Run: matlab -batch "make_sec4_figures"

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(this_dir, 'data'); fig_dir = fullfile(this_dir, 'figures'); if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
set(groot, 'defaultAxesFontName', 'Times New Roman', 'defaultTextFontName', 'Times New Roman', 'defaultAxesFontSize', 9, ...
    'defaultLegendFontSize', 7, 'defaultAxesTickDir', 'in', 'defaultAxesBox', 'on', 'defaultLineLineWidth', 0.9);
ylab_S = 'd^2W/d\omegad\Omega [J/(rad/s)/sr]';
N = 2000; xt = [15 20 30 40 50 70 100 150];

%% 4.1 undulator, cold bunch
U = readmatrix(fullfile(data_dir, 'undulator_N2000_onaxis.txt'), 'CommentStyle', '%');
f = figure('Units', 'centimeters', 'Position', [2 2 12 11]);
ax1 = subplot(2,1,1); hold on;
plot(U(:,1), U(:,3), 'k-', 'DisplayName', 'MAGNETO |\Sigma_jE_j|^2');
plot(U(:,1), U(:,5), 'r--', 'DisplayName', 'S_1[N + N(N-1)|F|^2]');
plot(U(:,1), N*U(:,4), ':', 'Color', [0 0.45 0.74], 'DisplayName', 'N S_1');
plot(U(:,1), N^2*U(:,4), ':', 'Color', [0.47 0.67 0.19], 'DisplayName', 'N^2 S_1');
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', xt); xlim([15 150]); ylabel(ylab_S); legend('Location', 'southeast');
text(0.02, 0.92, '(a)', 'Units', 'normalized');
ax2 = subplot(2,1,2); hold on;
plot(U(:,1), U(:,12), 'r-', 'DisplayName', 'analytic Gaussian form factor');
plot(U(:,1), U(:,13), '-', 'Color', [0.93 0.69 0.13], 'DisplayName', 'discrete form factor (same sample)');
plot(U(:,1), U(:,14), '-', 'Color', [0 0.45 0.74], 'DisplayName', 'trajectory translation');
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', xt); xlim([15 150]); ylim([1e-6 1]);
xlabel('\lambda [nm]'); ylabel('|\DeltaS| / envelope'); legend('Location', 'northeast');
text(0.02, 0.92, '(b)', 'Units', 'normalized'); linkaxes([ax1 ax2], 'x');
save_fig(f, fullfile(this_dir, 'figures_unused'), 'fig4_1_undulator_cold_two_panel');   % paper version: make_fig4_1.m (2026-10-28)

%% 4.2 dipole, cold bunch
D = readmatrix(fullfile(data_dir, 'dipole_N2000_onaxis.txt'), 'CommentStyle', '%');
SN = readmatrix(fullfile(data_dir, 'dipole_shot_noise.txt'), 'CommentStyle', '%');
f = figure('Units', 'centimeters', 'Position', [2 2 12 11]);
ax1 = subplot(2,1,1); hold on;
plot(D(:,1), D(:,3), 'k-', 'DisplayName', 'MAGNETO |\Sigma_jE_j|^2');
plot(D(:,1), D(:,5), 'r--', 'DisplayName', 'S_1[N + N(N-1)|F|^2]');
plot(D(:,1), N*D(:,4), ':', 'Color', [0 0.45 0.74], 'DisplayName', 'N S_1');
plot(D(:,1), N^2*D(:,4), ':', 'Color', [0.47 0.67 0.19], 'DisplayName', 'N^2 S_1');
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', xt); xlim([15 150]); ylabel(ylab_S); legend('Location', 'southeast');
text(0.02, 0.92, '(a)', 'Units', 'normalized');
ax2 = subplot(2,1,2); hold on;
fill([SN(:,1); flipud(SN(:,1))], [SN(:,2); -flipud(SN(:,2))], [0.85 0.85 0.85], 'EdgeColor', 'none', 'DisplayName', '\pm1\sigma shot noise (N = 2000)');
plot(D(:,1), D(:,17), 'k-', 'DisplayName', 'MAGNETO / analytic - 1');
yline(0, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
set(gca, 'XScale', 'log', 'XTick', xt, 'Layer', 'top'); xlim([15 150]); ylim([-1.2 1.2]);
xlabel('\lambda [nm]'); ylabel('relative deviation'); legend('Location', 'southeast');
text(0.02, 0.92, '(b)', 'Units', 'normalized'); linkaxes([ax1 ax2], 'x');
save_fig(f, fullfile(this_dir, 'figures_unused'), 'fig4_2_dipole_cold_two_panel');   % paper version: make_fig4_2.m (2026-10-28)

%% 4.3 - 4.5 energy spread / emittance
A = load(fullfile(data_dir, 'spread_emittance', 'sec4_spread_emittance.mat'));
C = A.cases_info; sd = [C.sig_d]; en = [C.en];
isd = find(en == 0);  [~, o] = sort(sd(isd)); isd = isd(o);
ien = find(sd == 0);  [~, o] = sort(en(ien)); ien = ien(o);
olab = {'on axis', '\theta_x = 0.5 mrad'};
fig_dir_main = fig_dir;
for q = 1:2
    on = A.oname{q}; hl = A.harm{q}; Wq = A.W.(on); Tq = A.T.(on); lam = A.lam;
    if q == 2; fig_dir = fullfile(this_dir, 'figures_unused'); else; fig_dir = fig_dir_main; end   % 0.5 mrad: not in the paper (2026-10-27)
    lam1 = Wq(1).lpk_inc(1);          % baseline n = 1 peak; lambda_n = lam1/n for the window centres

    % 4.3 incoherent line shapes
    f = figure('Units', 'centimeters', 'Position', [2 2 4.2*numel(hl) + 4 10]); tl = tiledlayout(2, numel(hl) + 1, 'TileSpacing', 'compact');
    for row = 1:2
        if row == 1; idx = isd; lab = arrayfun(@(k) sprintf('\\sigma_\\delta = %g%%', 100*sd(k)), idx, 'UniformOutput', false);
        else;        idx = ien; lab = arrayfun(@(k) sprintf('\\epsilon_n = %g \\mum', 1e6*en(k)), idx, 'UniformOutput', false); end
        cols = parula(numel(idx) + 1);
        for n = 1:numel(hl)
            nexttile((row-1)*(numel(hl) + 1) + n); hold on;
            lc = lam1/hl(n); sel = abs(lam/lc - 1) < 0.3;
            for j = 1:numel(idx)
                Sj = Wq(idx(j)).S_inc(:); pk = max(Sj(sel & abs(lam/lc - 1) < 0.1));
                plot(lam(sel)/lc - 1, Sj(sel)/pk, '-', 'Color', cols(j,:), 'LineWidth', 0.8);
            end
            xlim([-0.3 0.3]); ylim([0 1.15]);
            if row == 2; xlabel('\lambda/\lambda_n - 1'); end
            title(sprintf('n = %d', hl(n)), 'FontWeight', 'normal');
            if n == 1; ylabel('\Sigma|E_j|^2 / peak'); end
            if n == numel(hl); lg = legend(lab, 'FontSize', 6); lg.Layout.Tile = row*(numel(hl) + 1); end
        end
    end
    title(tl, sprintf('incoherent sum, %s (top: energy spread, bottom: emittance)', olab{q}), 'FontSize', 9);
    save_fig(f, fullfile(this_dir, 'figures_unused'), sprintf('fig4_3_lines_inc_%s', on));   % paper version: make_fig4_3_lines_n1.m (2026-10-28)

    % 4.4 FWHM
    f = figure('Units', 'centimeters', 'Position', [2 2 21 7.5]); cols = lines(numel(hl)); tiledlayout(1, 2, 'TileSpacing', 'compact');
    for p = 1:2
        nexttile; hold on;
        if p == 1; idx = isd; xv = 100*sd(idx); xl = '\sigma_\delta [%]'; else; idx = ien; xv = 1e6*en(idx); xl = '\epsilon_n [\mum]'; end
        xv(xv == 0) = min(xv(xv > 0))/3;                    % baseline plotted at the left edge
        for n = 1:numel(hl)
            si = Tq.show_inc(idx, n)'; sc = Tq.show_coh(idx, n)';
            fi = arrayfun(@(r) r.fwhm_inc(n), Wq(idx)); ei = arrayfun(@(r) r.err_inc(n), Wq(idx));
            fc = arrayfun(@(r) r.fwhm_coh(n), Wq(idx)); fe = arrayfun(@(r) r.fwhm_est(n), Wq(idx));
            fq = Tq.fwhm_quad(idx, n)';                       % textbook quadrature sum (on axis only, NaN otherwise)
            errorbar(xv(si), fi(si), ei(si), 'o-', 'Color', cols(n,:), 'MarkerFaceColor', cols(n,:), 'MarkerSize', 4, 'HandleVisibility', 'off');
            plot(xv(sc), fc(sc), 's', 'Color', cols(n,:), 'MarkerSize', 7, 'HandleVisibility', 'off');
            plot(xv(si), fe(si), '--', 'Color', cols(n,:), 'HandleVisibility', 'off');
            plot(xv(si), fq(si), ':', 'Color', cols(n,:), 'LineWidth', 1.4, 'HandleVisibility', 'off');
        end
        set(gca, 'XScale', 'log', 'YScale', 'log'); ylim([1e-2 0.6]); xlabel(xl); ylabel('FWHM / \lambda_{peak}');
        if p == 2                                             % colour = harmonic, style = quantity
            hh = gobjects(0); ll = {};
            for n = 1:numel(hl); hh(end+1) = plot(nan, nan, '-', 'Color', cols(n,:), 'LineWidth', 1.5); ll{end+1} = sprintf('n = %d', hl(n)); end %#ok<SAGROW>
            hh(end+1) = plot(nan, nan, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 4); ll{end+1} = '\Sigma|E_j|^2 (MAGNETO)';
            hh(end+1) = plot(nan, nan, 'ks', 'MarkerSize', 7);                         ll{end+1} = sprintf('|\\SigmaE_j|^2 (MAGNETO, r \\geq %g)', A.r_min);
            hh(end+1) = plot(nan, nan, 'k--');                                         ll{end+1} = 'simple estimate';
            if any(~isnan(Tq.fwhm_quad(:))); hh(end+1) = plot(nan, nan, 'k:', 'LineWidth', 1.4); ll{end+1} = 'quadrature formula'; end
            lg = legend(hh, ll, 'FontSize', 6); lg.Layout.Tile = 'east';
        end
        title(sprintf('%s (leftmost: none)', olab{q}), 'FontWeight', 'normal');
    end
    save_fig(f, fullfile(this_dir, 'figures_unused'), sprintf('fig4_4_fwhm_%s', on));   % not in the paper (user, 2026-10-28)

    % 4.5 coherence ratio
    f = figure('Units', 'centimeters', 'Position', [2 2 21 7.5]); tiledlayout(1, 2, 'TileSpacing', 'compact');
    for p = 1:2
        nexttile; hold on;
        if p == 1; idx = isd; xv = 100*sd(idx); xl = '\sigma_\delta [%]'; else; idx = ien; xv = 1e6*en(idx); xl = '\epsilon_n [\mum]'; end
        xv(xv == 0) = min(xv(xv > 0))/3;
        for n = 1:numel(hl)
            plot(xv, Tq.r_FF_d(idx, n), '-', 'Color', cols(n,:), 'HandleVisibility', 'off');
            plot(xv, Tq.r_FF_g(idx, n), '--', 'Color', cols(n,:), 'HandleVisibility', 'off');
            plot(xv, Tq.r_int(idx, n), 'o', 'Color', cols(n,:), 'MarkerFaceColor', cols(n,:), 'MarkerSize', 4, 'DisplayName', sprintf('n = %d, MAGNETO', hl(n)));
        end
        plot(nan, nan, 'k-', 'DisplayName', 'form factor, this sample');
        plot(nan, nan, 'k--', 'DisplayName', 'form factor, Gaussian \sigma_z = 10 nm');
        yline(A.r_min, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
        yline(1, '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
        set(gca, 'XScale', 'log', 'YScale', 'log'); ylim([5e-3 5e3]); xlabel(xl); ylabel('r = \int|\SigmaE_j|^2d\omega / \int\Sigma|E_j|^2d\omega (natural line)');
        if p == 2; lg = legend('FontSize', 6); lg.Layout.Tile = 'east'; end
        title(sprintf('%s (leftmost: none)', olab{q}), 'FontWeight', 'normal');
    end
    save_fig(f, fullfile(this_dir, 'figures_unused'), sprintf('fig4_5_coherence_%s', on));   % not in the paper (user, 2026-10-28)
end

function save_fig(f, d, name)
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(d, [name '.png']), 'Resolution', 300);
exportgraphics(f, fullfile(d, [name '.pdf']), 'ContentType', 'vector');
savefig(f, fullfile(d, [name '.fig'])); close(f);
end
