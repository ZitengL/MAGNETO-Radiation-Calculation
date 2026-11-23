%% plot_y2mm
% Spectra d2W/dw dOmega at the observer x = 0, y = -2.0 mm, z = 5.9 m (user's choice, 2026-11-05) for every sigma_z,
% from the 5e4 run (results/combined.mat via results/contrib.mat of contributions.m):
%   figure 1: no polariser (x+y+z): undulator part, dipole part, total (und + dip, incl. interference)
%   figure 2: x polariser: undulator x, dipole x, total x (passed) and the blocked light (y+z, total)
% A curve is drawn only where its jackknife error (25 chunks) is <= 30 %; exception: the dipole x part, which is
% below the macro-particle noise everywhere, is drawn as a thin dotted line labelled as noise level.
% Grey dotted: self terms sum_j |E_j|^2 of the total (expected macro-particle noise; the real incoherent level is 1/nmacro of it).
% Output: figures/spectra_y-2.0mm_unpol.png/.fig, figures/spectra_y-2.0mm_xpol.png/.fig, results/band_y-2.0mm.txt

clear; close all;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; eps0 = 8.85418781e-12; thr = 0.3; yo = -2.0;
load(fullfile(res_dir, 'contrib.mat'), 'S', 'lnm', 'sz', 'ys');
R = load(fullfile(res_dir, 'combined.mat')); R = R.R;
o = find(abs(ys - yo) < 1e-9); ns = numel(sz);
P = eps0*c0*norm(R.obs(o,:))^2/pi;
Iall = squeeze(sum(R.Ic, 1));                                       % obs x part x comp x omega
noise = {P*squeeze(sum(Iall(o,3,:,:), 3))', P*squeeze(Iall(o,3,1,:))'};  % self terms: unpolarised, x
g = @(p, f, k) squeeze(S(p).(f)(k,o,:))';
tk = [50 70 100 150 200 300];

for fig = 1:2
    f = figure('Position', [30 30 1500 1150], 'Visible', 'off');
    tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    if fig == 1
        title(tl, sprintf('x = 0, y = %.1f mm, z = 5.9 m, no polariser (x+y+z), 5e4 macro-particles; curves shown where the jackknife error <= 30 %%', yo));
    else
        title(tl, sprintf('x = 0, y = %.1f mm, z = 5.9 m, x polariser, 5e4 macro-particles; curves shown where the jackknife error <= 30 %%', yo));
    end
    allv = [];
    for k = 1:ns
        ax(k) = nexttile; hold on; %#ok<SAGROW>
        if fig == 1
            p = 1;
            c = {mk(g(p,'U',k), g(p,'sdU',k)), mk(g(p,'D',k), g(p,'sdD',k)), mk(g(p,'T',k), g(p,'sdT',k))};
            plot(lnm, c{1}, 'b-', 'LineWidth', 1.3, 'DisplayName', 'undulator part');
            plot(lnm, c{2}, 'r-', 'LineWidth', 1.3, 'DisplayName', 'dipole part');
            plot(lnm, c{3}, 'k-', 'LineWidth', 1.8, 'DisplayName', 'total (und + dip)');
        else
            c = {mk(g(2,'U',k), g(2,'sdU',k)), g(2,'D',k), mk(g(2,'T',k), g(2,'sdT',k)), mk(g(3,'T',k), g(3,'sdT',k))};
            plot(lnm, c{1}, 'b-', 'LineWidth', 1.3, 'DisplayName', 'undulator part, x');
            plot(lnm, c{2}, 'r:', 'LineWidth', 1.0, 'DisplayName', 'dipole part, x (macro-particle noise level)');
            plot(lnm, c{3}, 'k-', 'LineWidth', 1.8, 'DisplayName', 'total x (passed)');
            plot(lnm, c{4}, 'm--', 'LineWidth', 1.2, 'DisplayName', 'total y+z (blocked)');
        end
        plot(lnm, noise{fig}, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'DisplayName', 'self terms (noise)');
        set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', tk); grid on; box on; xlim([50 300]);
        xlabel('\lambda [nm]'); ylabel('d^2W/d\omegad\Omega [J/(rad/s)/sr]');
        title(sprintf('\\sigma_z = %g nm', sz(k)));
        if k == 1; legend('Location', 'southeast', 'FontSize', 7); end
        for i = 1:numel(c); v = c{i}; allv = [allv v(isfinite(v) & v > 0)]; end %#ok<AGROW>
        allv = [allv noise{fig}]; %#ok<AGROW>
    end
    yl = [10^floor(log10(max(min(allv), max(allv)*1e-8))) 10^ceil(log10(max(allv)))];
    set(ax, 'YLim', yl);
    nm = {'unpol', 'xpol'};
    saveas(f, fullfile(fig_dir, sprintf('spectra_y%.1fmm_%s.png', yo, nm{fig}))); savefig(f, fullfile(fig_dir, sprintf('spectra_y%.1fmm_%s.fig', yo, nm{fig})));
end

% band energies 50-300 nm
om = R.omega; Wn = @(v) -trapz(om, v);
fid = fopen(fullfile(res_dir, sprintf('band_y%.1fmm.txt', yo)), 'w');
fprintf(fid, '%% x = 0, y = %.1f mm, z = 5.9 m, energy in 50-300 nm [J/sr], 5e4 macro-particles; sd = jackknife (25 chunks)\n', yo);
fprintf(fid, '%% unpolarised: W_und W_dip W_tot W_tot_sd | x polariser: W_und_x W_dip_x(noise level) W_tot_x W_tot_x_sd W_self_x | blocked: W_tot_yz\n');
fprintf(fid, '%% sigma_z_nm  unpol_und unpol_dip unpol_tot unpol_tot_sd  x_und x_dip x_tot x_tot_sd x_self  yz_tot\n');
for k = 1:ns
    fprintf(fid, '%3.0f  %.4e %.4e %.4e %.2e  %.4e %.4e %.4e %.2e %.4e  %.4e\n', sz(k), S(1).WU(k,o), S(1).WD(k,o), S(1).WT(k,o), S(1).sdWT(k,o), ...
        S(2).WU(k,o), S(2).WD(k,o), S(2).WT(k,o), S(2).sdWT(k,o), Wn(noise{2}), S(3).WT(k,o));
end
fclose(fid); type(fullfile(res_dir, sprintf('band_y%.1fmm.txt', yo)));

function v = mk(v, sd)
% hide the points whose jackknife error is > 30 %
v(sd > 0.3*v) = nan;
end
