%% plot_y2mm_full
% "Complete" version of plot_y2mm.m at x = 0, y = -2.0 mm, z = 5.9 m (user's request, 2026-11-05):
%  - spectra over the whole 50-300 nm for every sigma_z, nothing hidden; shaded band = +- 1 jackknife standard
%    deviation (25 chunks); grey dotted = self terms (macro-particle noise level);
%  - energy in 50-300 nm vs sigma_z. The band energy W integrates the WHOLE 50-300 nm, i.e. it includes the
%    noise-dominated short wavelengths; to show how much they matter, W_x is also given
%      (a) only over the resolved wavelengths (jackknife error <= 30 %),
%      (b) with the self terms subtracted (removes the mean noise bias, not the fluctuation).
% Output: figures/spectra_y-2.0mm_unpol_full.png/.fig, figures/spectra_y-2.0mm_xpol_full.png/.fig,
%         figures/W_vs_sigz_y-2.0mm.png/.fig, results/W_vs_sigz_y-2.0mm.txt

clear; close all;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; eps0 = 8.85418781e-12; thr = 0.3; yo = -2.0;
load(fullfile(res_dir, 'contrib.mat'), 'S', 'lnm', 'sz', 'ys');
R = load(fullfile(res_dir, 'combined.mat')); R = R.R; om = R.omega;
o = find(abs(ys - yo) < 1e-9); ns = numel(sz);
P = eps0*c0*norm(R.obs(o,:))^2/pi;
Iall = squeeze(sum(R.Ic, 1));                                        % obs x part x comp x omega
noise = {P*squeeze(sum(Iall(o,3,:,:), 3))', P*squeeze(Iall(o,3,1,:))', P*squeeze(sum(Iall(o,3,2:3,:), 3))'};  % unpol, x, y+z
g = @(p, f, k) squeeze(S(p).(f)(k,o,:))';
Wn = @(v) -trapz(om, v);
tk = [50 70 100 150 200 300];

%% spectra, nothing hidden
for fig = 1:2
    f = figure('Position', [30 30 1500 1150], 'Visible', 'off');
    tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    if fig == 1; s1 = 'no polariser (x+y+z)'; else; s1 = 'x polariser'; end
    title(tl, sprintf('x = 0, y = %.1f mm, z = 5.9 m, %s, 5e4 macro-particles, complete 50-300 nm; band = \\pm1 jackknife s.d.', yo, s1));
    allv = [];
    for k = 1:ns
        ax(k) = nexttile; hold on; %#ok<SAGROW>
        if fig == 1
            cur = {1,'U','b','-',1.3,'undulator part'; 1,'D','r','-',1.3,'dipole part'; 1,'T','k','-',1.8,'total (und + dip)'};
        else
            cur = {2,'U','b','-',1.3,'undulator part, x'; 2,'D','r',':',1.0,'dipole part, x'; 2,'T','k','-',1.8,'total x (passed)'; 3,'T','m','--',1.2,'total y+z (blocked)'};
        end
        for i = 1:size(cur, 1)
            v = g(cur{i,1}, cur{i,2}, k); sd = g(cur{i,1}, ['sd' cur{i,2}], k);
            lo = max(v - sd, v/1e3); hi = v + sd;
            fill([lnm fliplr(lnm)], [lo fliplr(hi)], cur{i,3}, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(lnm, v, 'Color', cur{i,3}, 'LineStyle', cur{i,4}, 'LineWidth', cur{i,5}, 'DisplayName', cur{i,6});
            allv = [allv v]; %#ok<AGROW>
        end
        plot(lnm, noise{fig}, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2, 'DisplayName', 'self terms (noise level)');
        allv = [allv noise{fig}]; %#ok<AGROW>
        set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', tk, 'Layer', 'top'); grid on; box on; xlim([50 300]);
        xlabel('\lambda [nm]'); ylabel('d^2W/d\omegad\Omega [J/(rad/s)/sr]'); title(sprintf('\\sigma_z = %g nm', sz(k)));
        if k == 1; legend('Location', 'southeast', 'FontSize', 7); end
    end
    yl = [10^floor(log10(max(min(allv(allv > 0)), max(allv)*1e-8))) 10^ceil(log10(max(allv)))];
    set(ax, 'YLim', yl);
    nm = {'unpol', 'xpol'};
    saveas(f, fullfile(fig_dir, sprintf('spectra_y%.1fmm_%s_full.png', yo, nm{fig}))); savefig(f, fullfile(fig_dir, sprintf('spectra_y%.1fmm_%s_full.fig', yo, nm{fig})));
end

%% band energy vs sigma_z
Wself = cellfun(Wn, noise);                                          % unpol, x, y+z
T = zeros(ns, 12);
for k = 1:ns
    v = g(2,'T',k); sd = g(2,'sdT',k); ok = sd <= thr*v;
    lmin = min(lnm(ok));
    T(k,:) = [sz(k), S(1).WU(k,o), S(1).WD(k,o), S(1).WT(k,o), S(1).sdWT(k,o), ...
              S(2).WU(k,o), S(2).WT(k,o), S(2).sdWT(k,o), Wn(v.*ok), S(2).WT(k,o) - Wself(2), lmin, S(3).WT(k,o)];
end
fid = fopen(fullfile(res_dir, sprintf('W_vs_sigz_y%.1fmm.txt', yo)), 'w');
fprintf(fid, '%% x = 0, y = %.1f mm, z = 5.9 m, energy [J/sr]; W = integral over the WHOLE 50-300 nm (incl. noise-dominated wavelengths); sd = jackknife\n', yo);
fprintf(fid, '%% self terms (macro-particle noise level, independent of sigma_z): unpol %.4e, x %.4e, y+z %.4e\n', Wself);
fprintf(fid, '%% W_x_resolved = integral only where the jackknife error of S_tot,x <= 30 %% (lambda >= lambda_min_x); W_x_minus_self = W_tot_x - self_x\n');
fprintf(fid, '%% sigma_z_nm unpol_W_und unpol_W_dip unpol_W_tot unpol_W_tot_sd x_W_und x_W_tot x_W_tot_sd x_W_resolved x_W_minus_self lambda_min_x_nm yz_W_tot\n');
fprintf(fid, '%3.0f %.4e %.4e %.4e %.2e %.4e %.4e %.2e %.4e %.4e %6.1f %.4e\n', T');
fclose(fid); type(fullfile(res_dir, sprintf('W_vs_sigz_y%.1fmm.txt', yo)));

f = figure('Position', [30 30 1400 560], 'Visible', 'off');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('x = 0, y = %.1f mm, z = 5.9 m: energy in 50-300 nm (whole band) vs \\sigma_z, 5e4 macro-particles, error bars = jackknife', yo));
nexttile; hold on;
plot(sz, T(:,2), 'b-o', 'DisplayName', 'undulator part');
plot(sz, T(:,3), 'r-s', 'DisplayName', 'dipole part');
errorbar(sz, T(:,4), T(:,5), 'k-o', 'LineWidth', 1.5, 'DisplayName', 'total (und + dip)');
yline(Wself(1), ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2, 'DisplayName', 'self terms (noise level)');
set(gca, 'YScale', 'log'); grid on; box on; xlabel('\sigma_z [nm]'); ylabel('energy in 50-300 nm [J/sr]'); title('no polariser (x+y+z)');
legend('Location', 'southwest');
nexttile; hold on;
errorbar(sz, T(:,7), T(:,8), 'k-o', 'LineWidth', 1.5, 'DisplayName', 'total x (passed), whole band');
plot(sz, T(:,9), 'g--^', 'LineWidth', 1.2, 'DisplayName', 'total x, resolved \lambda only');
plot(sz, T(:,10), 'c-.d', 'LineWidth', 1.2, 'DisplayName', 'total x minus self terms');
plot(sz, T(:,6), 'b:o', 'DisplayName', 'undulator part, x');
plot(sz, T(:,12), 'm--s', 'DisplayName', 'total y+z (blocked)');
yline(Wself(2), ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2, 'DisplayName', 'self terms x (noise level)');
set(gca, 'YScale', 'log'); grid on; box on; xlabel('\sigma_z [nm]'); ylabel('energy in 50-300 nm [J/sr]'); title('x polariser');
legend('Location', 'southwest');
saveas(f, fullfile(fig_dir, sprintf('W_vs_sigz_y%.1fmm.png', yo))); savefig(f, fullfile(fig_dir, sprintf('W_vs_sigz_y%.1fmm.fig', yo)));
