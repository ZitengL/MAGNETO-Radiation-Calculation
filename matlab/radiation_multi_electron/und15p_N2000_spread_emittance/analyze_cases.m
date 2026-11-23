%% analyze_cases
% Energy-spread / emittance study (run_all.m): broadening of the undulator harmonics.
% For every case and observer: coherent spectrum S_coh = |sum E_j|^2 and incoherent sum S_inc = sum |E_j|^2
% (pref = eps0 c R^2/pi, per sr of the 11 m sphere). For each harmonic:
%   peak position (parabola through the 3 highest points, in log lambda) and FWHM of the main lobe
%   (half-maximum crossings on both sides, linear interpolation), relative to the peak wavelength;
%   statistical error: the same from the two halves (odd / even IDs, 1000 electrons each),
%   error of the full value ~ |FWHM_odd - FWHM_even|/2.
% Simple estimate (reference only, not a theory): every electron radiates the baseline line shape,
% shifted to its own resonance: d lambda/lambda = -2 delta_j + gamma^2 (theta_j^2 - theta_o^2)/(1 + K^2/2 + gamma^2 theta_o^2),
% theta_j = angle between the electron direction and the observer direction seen from the undulator
% centre, theta_o = that angle for an ideal electron; S_est(lambda) = mean_j S_inc,baseline(lambda/(1 + eps_j)).
% Outputs: results/analysis.mat, results/widths.txt, figures/*.png/.fig

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E240MeV_und15p_N2000_spread_emittance'); addpath(gdir);
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures'); if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
c0 = 299792458; eps0 = 8.85418781e-12; gam0 = 470.66442480; K = 0.980309; zu = 1.0;
C = cases(); C = C(arrayfun(@(c) exist(fullfile(res_dir, [c.name '.mat']), 'file') == 2, C));
fprintf('%d cases available\n', numel(C));
R0 = load(fullfile(res_dir, 'sd0_en0.mat')); lam = R0.res.lambda_nm; obs = R0.res.obs; oname = R0.res.obs_name;
olab = {'on axis', '\theta_x = 0.5 mrad'};
harm = {[1 3 5], [1 2 4 5]};
lam1 = [101.8283, 101.8283*(1 + gam0^2*atan(11*sin(0.5e-3)/(11*cos(0.5e-3) - 1))^2/(1 + K^2/2))];   % lambda_1(theta) [nm]
spec = @(E, q, h, pref) pref*squeeze(sum(abs(E(q, h, :, :)).^2, 3))';
spinc = @(I, q, h, pref) pref*squeeze(sum(I(q, h, :, :), 3))';

W = struct();
for q = 1:2
    pref = eps0*c0*norm(obs{q})^2/pi;
    th_o = atan(obs{q}(1)/(obs{q}(3) - zu));
    Sb = spinc(R0.res.I, q, 1, pref);
    for i = 1:numel(C)
        R = load(fullfile(res_dir, [C(i).name '.mat'])); R = R.res;
        D = readmatrix(fullfile(gdir, C(i).name, 'initial_distribution.txt'), 'NumHeaderLines', 1);
        xp = D(:,4)./D(:,6); yp = D(:,5)./D(:,6); g = D(:,7); dl = g/gam0 - 1;
        xw = D(:,1) + xp*zu; yw = D(:,2) + yp*zu;                               % position at the waist (undulator centre)
        ox = atan((obs{q}(1) - xw)/(obs{q}(3) - zu)); oy = atan((obs{q}(2) - yw)/(obs{q}(3) - zu));
        thj2 = (ox - xp).^2 + (oy - yp).^2;
        epsj = -2*dl + gam0^2*(thj2 - th_o^2)/(1 + K^2/2 + gam0^2*th_o^2);
        Sest = zeros(size(lam));
        for j = 1:numel(epsj); Sest = Sest + interp1(lam, Sb, lam/(1 + epsj(j)), 'linear', 0); end
        S = {spec(R.E, q, 1, pref), spinc(R.I, q, 1, pref), Sest, spec(R.E, q, 2, pref), spec(R.E, q, 3, pref), spinc(R.I, q, 2, pref), spinc(R.I, q, 3, pref)};
        hl = harm{q}; r = struct('n', hl);
        for n = 1:numel(hl)
            lc = lam1(q)/hl(n);
            for m = 1:numel(S)
                [lp(m), fw(m), pk(m)] = peak_fwhm(lam, S{m}, lc); %#ok<SAGROW>
            end
            r.lpk_coh(n) = lp(1); r.fwhm_coh(n) = fw(1); r.pk_coh(n) = pk(1);
            r.lpk_inc(n) = lp(2); r.fwhm_inc(n) = fw(2); r.pk_inc(n) = pk(2);
            r.lpk_est(n) = lp(3); r.fwhm_est(n) = fw(3);
            r.err_coh(n) = abs(fw(4) - fw(5))/2; r.err_inc(n) = abs(fw(6) - fw(7))/2;
        end
        r.S_coh = S{1}; r.S_inc = S{2}; r.S_est = Sest; r.sig_d = C(i).sig_d; r.en = C(i).en; r.name = C(i).name;
        r.sample_sig_d = std(dl); r.sample_sig_th = std(xp); r.xmax = R.xmax_in_map;
        W.(oname{q})(i) = r;
    end
end

%% table
fid = fopen(fullfile(res_dir, 'widths.txt'), 'wt');
fprintf(fid, '%% Energy-spread / emittance study, 2000 electrons, 15-period undulator, 240 MeV. FWHM relative to the peak wavelength.\n');
fprintf(fid, '%% inc = incoherent sum sum|E_j|^2, coh = coherent |sum E_j|^2, est = simple estimate (resonance shift of each electron); err = half the difference of the two halves\n');
for q = 1:2
    fprintf(fid, '\n%% observer %s\n%% case            sig_d    en[um]  n  peak_inc[nm]  FWHM_inc   err_inc    FWHM_coh   err_coh    FWHM_est   S_coh/S_inc(peak)  max|x,y| in map [mm]\n', oname{q});
    fprintf('\n=== %s\n', oname{q});
    for i = 1:numel(C)
        r = W.(oname{q})(i);
        for n = 1:numel(r.n)
            line = sprintf('%-14s  %-7.4g  %-6.3g  %d  %11.4f  %9.3e  %9.2e  %9.3e  %9.2e  %9.3e  %12.3e  %8.3f', r.name, r.sig_d, 1e6*r.en, r.n(n), r.lpk_inc(n), r.fwhm_inc(n), r.err_inc(n), ...
                r.fwhm_coh(n), r.err_coh(n), r.fwhm_est(n), r.pk_coh(n)/r.pk_inc(n), 1e3*r.xmax);
            fprintf(fid, '%s\n', line); fprintf('%s\n', line);
        end
    end
end
fclose(fid);
save(fullfile(res_dir, 'analysis.mat'), 'W', 'lam', 'harm', 'lam1');

%% figures
isd = find([C.en] == 0); [~, o] = sort([C(isd).sig_d]); isd = isd(o);
ien = find([C.sig_d] == 0); [~, o] = sort([C(ien).en]); ien = ien(o);
for q = 1:2
    hl = harm{q}; Wq = W.(oname{q});
    for kind = 1:2
        f = figure('Position', [50 50 300*numel(hl) + 200 750]);
        for row = 1:2
            if row == 1; idx = isd; lab = arrayfun(@(c) sprintf('\\sigma_\\delta = %g%%', 100*c.sig_d), C(idx), 'UniformOutput', false);
            else;        idx = ien; lab = arrayfun(@(c) sprintf('\\epsilon_n = %g \\mum', 1e6*c.en), C(idx), 'UniformOutput', false); end
            cols = parula(numel(idx) + 1);
            for n = 1:numel(hl)
                subplot(2, numel(hl), (row-1)*numel(hl) + n); hold on;
                lc = lam1(q)/hl(n); sel = abs(lam/lc - 1) < 0.3;
                for j = 1:numel(idx)
                    if kind == 1; Sj = Wq(idx(j)).S_inc; else; Sj = Wq(idx(j)).S_coh; end
                    pk = max(Sj(sel & abs(lam/lc - 1) < 0.1));
                    plot(lam(sel)/lc - 1, Sj(sel)/pk, '-', 'Color', cols(j,:), 'LineWidth', 1);
                end
                grid on; xlim([-0.3 0.3]); ylim([0 1.15]);
                title(sprintf('n = %d (\\lambda_n = %.2f nm)', hl(n), lc)); xlabel('\lambda/\lambda_n - 1');
                if n == 1; ylabel('normalised to the peak'); legend(lab, 'Location', 'northwest', 'FontSize', 7); end
            end
        end
        if kind == 1; kt = 'incoherent sum \Sigma|E_j|^2'; tg = 'inc'; else; kt = 'coherent |\SigmaE_j|^2'; tg = 'coh'; end
        sgtitle(sprintf('2000 electrons, 15-period undulator, 240 MeV, %s: %s (top: energy spread, bottom: emittance)', olab{q}, kt));
        set(findall(f, 'Type', 'axes'), 'Toolbar', []);
        exportgraphics(f, fullfile(fig_dir, sprintf('lines_%s_%s.png', tg, oname{q})), 'Resolution', 200); savefig(f, fullfile(fig_dir, sprintf('lines_%s_%s.fig', tg, oname{q}))); close(f);
    end
    % FWHM vs parameter
    f = figure('Position', [50 50 1100 450]); cols = lines(numel(hl));
    for p = 1:2
        subplot(1,2,p); hold on;
        if p == 1; idx = isd; xv = 100*[C(idx).sig_d]; xl = '\sigma_\delta [%]'; else; idx = ien; xv = 1e6*[C(idx).en]; xl = '\epsilon_n [\mum]'; end
        if any(xv > 0); xv(xv == 0) = min(xv(xv > 0))/3; else; xv(xv == 0) = 1e-2; end   % baseline plotted at the left edge
        for n = 1:numel(hl)
            fi = arrayfun(@(r) r.fwhm_inc(n), Wq(idx)); ei = arrayfun(@(r) r.err_inc(n), Wq(idx));
            fc = arrayfun(@(r) r.fwhm_coh(n), Wq(idx)); fe = arrayfun(@(r) r.fwhm_est(n), Wq(idx));
            errorbar(xv, fi, ei, 'o-', 'Color', cols(n,:), 'MarkerFaceColor', cols(n,:), 'DisplayName', sprintf('n = %d, incoherent', hl(n)));
            plot(xv, fc, 's', 'Color', cols(n,:), 'DisplayName', sprintf('n = %d, coherent', hl(n)));
            plot(xv, fe, '--', 'Color', cols(n,:), 'DisplayName', sprintf('n = %d, simple estimate', hl(n)));
        end
        set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; xlabel(xl); ylabel('FWHM / \lambda_{peak}');
        if p == 1; legend('Location', 'northwest', 'FontSize', 7); end
        title(sprintf('%s (leftmost point: no spread / emittance)', olab{q}));
    end
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(f, fullfile(fig_dir, sprintf('fwhm_%s.png', oname{q})), 'Resolution', 200); savefig(f, fullfile(fig_dir, sprintf('fwhm_%s.fig', oname{q}))); close(f);
    % overview spectra
    f = figure('Position', [50 50 1000 700]); sel = {'sd0_en0', 'sd3_en0', 'sd0_en12', 'sd1_en3'}; cl = lines(4);
    for kind = 1:2
        subplot(2,1,kind); hold on;
        for j = 1:numel(sel)
            i = find(strcmp({C.name}, sel{j})); if isempty(i); continue; end
            if kind == 1; Sj = Wq(i).S_coh; else; Sj = Wq(i).S_inc; end
            plot(lam, Sj, '-', 'Color', cl(j,:), 'DisplayName', strrep(sel{j}, '_', ', '));
        end
        set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; xlim([15 150]); set(gca, 'XTick', [15 20 30 40 50 70 100 150]);
        ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]'); legend('Location', 'southeast');
        if kind == 1; title(['coherent |\SigmaE_j|^2, ' olab{q}]); else; title(['incoherent \Sigma|E_j|^2, ' olab{q}]); xlabel('\lambda [nm]'); end
    end
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(f, fullfile(fig_dir, sprintf('overview_%s.png', oname{q})), 'Resolution', 200); savefig(f, fullfile(fig_dir, sprintf('overview_%s.fig', oname{q}))); close(f);
end

function [lp, fw, pk] = peak_fwhm(lam, S, lc)
% peak near lc (within +-10%), refined; FWHM of the main lobe (contiguous above half maximum), relative
win = find(abs(lam/lc - 1) < 0.1); [~, j] = max(S(win)); j = win(j);
x = log(lam(j-1:j+1)); y = S(j-1:j+1); p = polyfit(x - x(2), y, 2); dx = -p(2)/(2*p(1));
if abs(dx) > abs(x(3) - x(2)); dx = 0; end
lp = exp(x(2) + dx); pk = polyval(p, dx); h = pk/2;
a = j; while a > 1 && S(a) > h; a = a - 1; end
b = j; while b < numel(S) && S(b) > h; b = b + 1; end
la = interp1(S([a a+1]), lam([a a+1]), h); lb = interp1(S([b-1 b]), lam([b-1 b]), h);
fw = (lb - la)/lp;
end
