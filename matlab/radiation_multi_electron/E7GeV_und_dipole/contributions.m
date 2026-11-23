%% contributions
% Contributions of the undulator part and of the dipole part to the total radiation, for every sigma_z and wavelength,
% from the 5e4 run (results/combined.mat), observers x = 0, y = -1 ... -3 mm, z = 5.9 m.
%   S_und = pref |sum E_und|^2, S_dip = pref |sum E_dip|^2, S_tot = pref |sum (E_und + E_dip)|^2,
%   interference = S_tot - S_und - S_dip; summed over the field components of each polarisation:
%   unpolarised = x + y + z, passed by the x polariser = x, blocked = y + z.
% Error: jackknife over the 25 chunks (leave one chunk out, fields x 25/24). A value is shown only when its
%   jackknife standard deviation is <= 30 % of its magnitude (thr); otherwise it is grey (= not resolved by 5e4
%   macro-particles; for the dipole x part this is the macro-particle noise).
% Output: figures/contrib_unpolarised.png/.fig, figures/contrib_xpol.png/.fig, figures/polariser_passed_vs_blocked.png/.fig,
%         results/contrib_band.txt, results/contrib_spectra_y*.txt, results/contrib.mat

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; eps0 = 8.85418781e-12; thr = 0.3;
R = load(fullfile(res_dir, 'combined.mat')); R = R.R;
om = R.omega; lnm = R.lambda*1e9; sz = R.sigma_z*1e9; ns = numel(sz); no = size(R.obs, 1); nw = numel(om);
ys = 1e3*R.obs(:,2)'; nch = size(R.Ec, 1);
P = reshape(eps0*c0*vecnorm(R.obs, 2, 2).^2/pi, [1 1 no 1]);
Ea = sum(R.Ec, 1);                                  % 1 x sz x obs x part x comp x omega
Ej = (Ea - R.Ec)*nch/(nch-1);                       % leave-one-chunk-out
jsd = @(X) reshape(sqrt((nch-1)/nch*sum((X - mean(X, 1)).^2, 1)), [ns no size(X, 4)]);
sq = @(X) reshape(X, [ns no size(X, 4)]);
Wb = @(X) -trapz(om, X, 4);                         % band 50-300 nm

pol = struct('name', {'unpolarised (x+y+z)', 'x (passed by the polariser)', 'y+z (blocked)'}, 'tag', {'unpol', 'x', 'yz'}, 'comp', {1:3, 1, 2:3});
for p = 1:3
    [U, D, T] = parts(Ea, pol(p).comp, P); [Uj, Dj, Tj] = parts(Ej, pol(p).comp, P);
    s = struct();
    s.U = sq(U); s.D = sq(D); s.T = sq(T);
    s.sdU = jsd(Uj); s.sdD = jsd(Dj); s.sdT = jsd(Tj);
    s.rU = sq(U./T); s.rD = sq(D./T); s.rI = sq((T-U-D)./T);
    s.sdrU = jsd(Uj./Tj); s.sdrD = jsd(Dj./Tj); s.sdrI = jsd((Tj-Uj-Dj)./Tj);
    s.Q = sq(log10(T./(U+D))); s.sdQ = jsd(log10(Tj./(Uj+Dj)));       % interference: S_tot / (S_und + S_dip)
    % band energies
    s.WU = sq(Wb(U)); s.WD = sq(Wb(D)); s.WT = sq(Wb(T));
    s.sdWT = jsd(Wb(Tj));
    s.WrU = s.WU./s.WT; s.WrD = s.WD./s.WT; s.WrI = (s.WT - s.WU - s.WD)./s.WT;
    s.sdWrU = jsd(Wb(Uj)./Wb(Tj)); s.sdWrD = jsd(Wb(Dj)./Wb(Tj)); s.sdWrI = jsd((Wb(Tj) - Wb(Uj) - Wb(Dj))./Wb(Tj));
    S(p) = s; %#ok<SAGROW>
end

%% heat maps of the contributions (unpolarised and x)
bwr = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
tk = [50 70 100 150 200 300];
for p = 1:2
    s = S(p);
    f = figure('Position', [30 30 1500 1650], 'Visible', 'off');
    tl = tiledlayout(no, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, {sprintf('%s, x = 0, z = 5.9 m, 5e4 macro-particles', pol(p).name), ...
        sprintf('grey: jackknife error > %d %% (not resolved, macro-particle noise)', 100*thr)});
    C = {log10(s.rU), log10(s.rD), s.Q};
    OK = {s.sdrU <= thr*s.rU, s.sdrD <= thr*s.rD, s.sdQ <= log10(1 + thr)};
    lab = {'log_{10}(S_{und} / S_{tot})', 'log_{10}(S_{dip} / S_{tot})', 'log_{10}(S_{tot} / (S_{und} + S_{dip})): interference'};
    for k = 1:3
        v = C{k}(OK{k});
        if isempty(v); lim{k} = [-1 1];                                           %#ok<SAGROW>
        elseif k < 3; lim{k} = [floor(min(v)) ceil(max(v))];                       %#ok<SAGROW>
        else; a = ceil(2*max(abs(v)))/2; lim{k} = [-a a]; end                       %#ok<SAGROW>
        if diff(lim{k}) == 0; lim{k} = lim{k} + [-1 1]; end                        %#ok<SAGROW>
    end
    for o = 1:no
        for k = 1:3
            ax = nexttile;
            if k == 3; cm = bwr; else; cm = parula(256); end
            show_map(ax, squeeze(C{k}(:,o,:)), squeeze(OK{k}(:,o,:)), cm, lim{k}, log10(lnm), sz);
            set(ax, 'XTick', log10(tk), 'XTickLabel', string(tk));
            xlabel('\lambda [nm]'); ylabel('\sigma_z [nm]');
            title(sprintf('y = %.1f mm: %s', ys(o), lab{k}), 'FontSize', 9);
        end
    end
    saveas(f, fullfile(fig_dir, ['contrib_' pol(p).tag '.png'])); savefig(f, fullfile(fig_dir, ['contrib_' pol(p).tag '.fig']));
    for k = 1:3; fprintf('%s, column %d: colour range %g..%g, resolved cells %.1f %%\n', pol(p).tag, k, lim{k}, 100*mean(OK{k}(:))); end
end

%% passed vs blocked
f = figure('Position', [30 30 1600 900], 'Visible', 'off');
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'x polariser at x = 0, z = 5.9 m: passed (x, solid) and blocked (y+z, dashed); x shown only where the jackknife error <= 30 %');
sel = [1 3 5 9]; cl = lines(numel(sel));
for o = 1:no
    nexttile; hold on;
    for i = 1:numel(sel)
        X = squeeze(S(2).T(sel(i),o,:)); bad = squeeze(S(2).sdT(sel(i),o,:)) > thr*X; X(bad) = nan;
        Y = squeeze(S(3).T(sel(i),o,:)); badY = squeeze(S(3).sdT(sel(i),o,:)) > thr*Y; Y(badY) = nan;
        plot(lnm, X, '-', 'Color', cl(i,:), 'LineWidth', 1.3, 'DisplayName', sprintf('passed, \\sigma_z = %g nm', sz(sel(i))));
        plot(lnm, Y, '--', 'Color', cl(i,:), 'LineWidth', 1.0, 'DisplayName', sprintf('blocked, \\sigma_z = %g nm', sz(sel(i))));
    end
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', tk); grid on; box on; xlim([50 300]);
    xlabel('\lambda [nm]'); ylabel('d^2W/d\omegad\Omega [J/(rad/s)/sr]'); title(sprintf('y = %.1f mm', ys(o)));
    if o == 1; legend('Location', 'southwest', 'FontSize', 7); end
end
nexttile; hold on; co = lines(no);
for o = 1:no
    errorbar(sz, S(2).WT(:,o), S(2).sdWT(:,o), '-o', 'Color', co(o,:), 'DisplayName', sprintf('passed, y = %.1f mm', ys(o)));
    plot(sz, S(3).WT(:,o), '--s', 'Color', co(o,:), 'DisplayName', sprintf('blocked, y = %.1f mm', ys(o)));
end
set(gca, 'YScale', 'log'); grid on; box on; xlabel('\sigma_z [nm]'); ylabel('energy in 50-300 nm [J/sr]');
title('band energy'); legend('Location', 'southwest', 'FontSize', 7, 'NumColumns', 2);
saveas(f, fullfile(fig_dir, 'polariser_passed_vs_blocked.png')); savefig(f, fullfile(fig_dir, 'polariser_passed_vs_blocked.fig'));

%% tables
fid = fopen(fullfile(res_dir, 'contrib_band.txt'), 'w');
fprintf(fid, '%% band 50-300 nm [J/sr], 5e4 macro-particles, x = 0, z = 5.9 m; sd = jackknife over 25 chunks\n');
fprintf(fid, '%% pol: unpol = x+y+z, x = passed by the x polariser, yz = blocked\n');
fprintf(fid, '%% pol y_mm sigma_z_nm W_und W_dip W_tot W_tot_sd und/tot und/tot_sd dip/tot dip/tot_sd interference/tot interference/tot_sd\n');
for p = 1:3
    s = S(p);
    for o = 1:no
        for k = 1:ns
            fprintf(fid, '%-5s %5.1f %3.0f %.4e %.4e %.4e %.2e %.4e %.2e %.4e %.2e %+.4e %.2e\n', pol(p).tag, ys(o), sz(k), ...
                s.WU(k,o), s.WD(k,o), s.WT(k,o), s.sdWT(k,o), s.WrU(k,o), s.sdWrU(k,o), s.WrD(k,o), s.sdWrD(k,o), s.WrI(k,o), s.sdWrI(k,o));
        end
    end
end
fclose(fid);
for o = 1:no
    fid = fopen(fullfile(res_dir, sprintf('contrib_spectra_y%.1fmm.txt', ys(o))), 'w');
    fprintf(fid, '%% d2W/dw dOmega [J/(rad/s)/sr], y = %.1f mm; for each sigma_z (nm) and pol (unpol, x, yz): S_und S_dip S_tot S_tot_jackknife_sd\n%% lambda_nm', ys(o));
    for k = 1:ns; for p = 1:3; fprintf(fid, ' %s_und_%g %s_dip_%g %s_tot_%g %s_tot_sd_%g', pol(p).tag, sz(k), pol(p).tag, sz(k), pol(p).tag, sz(k), pol(p).tag, sz(k)); end; end
    fprintf(fid, '\n');
    for w = 1:nw
        fprintf(fid, '%.4f', lnm(w));
        for k = 1:ns; for p = 1:3; fprintf(fid, ' %.4e %.4e %.4e %.2e', S(p).U(k,o,w), S(p).D(k,o,w), S(p).T(k,o,w), S(p).sdT(k,o,w)); end; end
        fprintf(fid, '\n');
    end
    fclose(fid);
end
save(fullfile(res_dir, 'contrib.mat'), 'S', 'pol', 'lnm', 'sz', 'ys', 'thr');

%% numbers for the error discussion (x polarisation)
Iall = squeeze(sum(R.Ic, 1));                       % obs x part x comp x omega
for o = 1:no
    selfx = squeeze(P(o))*squeeze(Iall(o,3,1,:))'; Wself = -trapz(om, selfx);
    X = squeeze(S(2).T(ns,o,:))'; long = lnm >= 150;
    fprintf('y = %.1f mm, sigma_z = 50 nm, x: W = %.3e, self = %.3e (%.1f %%), share of W from lambda >= 150 nm = %.1f %%, jackknife = %.1f %%\n', ...
        ys(o), S(2).WT(ns,o), Wself, 100*Wself/S(2).WT(ns,o), 100*(-trapz(om(long), X(long)))/S(2).WT(ns,o), 100*S(2).sdWT(ns,o)/S(2).WT(ns,o));
end
Ez = squeeze(sum(abs(Ea(1,:,:,1,3,:) + Ea(1,:,:,2,3,:)).^2, 'all')); Ey = squeeze(sum(abs(Ea(1,:,:,1,2,:) + Ea(1,:,:,2,2,:)).^2, 'all'));
fprintf('blocked light: |E_z|^2 / |E_y|^2 summed over everything = %.2e\n', Ez/Ey);

% interference in x at a few wavelengths: 5e4 run (sigma_z = 10 nm) vs a single ideal electron
root_dir = fileparts(fileparts(fileparts(this_dir))); addpath(fullfile(root_dir, 'magneto'));
tr = read_gpt_trajectory(fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole', 'ideal_electron', 'result_sm_dz60.gdf'));
seg = {sub(tr, tr.r(3,:) < 0.63), sub(tr, tr.r(3,:) >= 0.63)};          % same split as run_combined.m
for o = 1:no
    Eu = emission_time_spectrum(seg{1}, R.obs(o,:), R.lambda, 'fB', norm(R.obs(o,:))/c0).Ew(1,:);
    Ed = emission_time_spectrum(seg{2}, R.obs(o,:), R.lambda, 'fB', norm(R.obs(o,:))/c0).Ew(1,:);
    Qi = abs(Eu+Ed).^2./(abs(Eu).^2 + abs(Ed).^2);
    for L = [100 200 300]
        [~, w] = min(abs(lnm - L));
        fprintf('y = %.1f mm, lambda %3d nm, x: S_tot/(S_und+S_dip) = %.3f +- %.3f (5e4, sigma_z = 10 nm), single ideal electron %.4f\n', ...
            ys(o), L, 10^S(2).Q(1,o,w), log(10)*10^S(2).Q(1,o,w)*S(2).sdQ(1,o,w), Qi(w));
    end
end

function t = sub(t, m)
k = find(m); kf = k(any(t.fB(:, k) ~= 0, 1)); k = max(k(1), kf(1)-1) : min(k(end), kf(end)+1);
for f = {'t', 'G'};             t.(f{1}) = t.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; t.(f{1}) = t.(f{1})(:, k); end
end

function show_map(ax, X, ok, cmap, lim, xv, yv)
% colour image with grey cells where ~ok (independent of the renderer's transparency)
n = size(cmap, 1);
idx = round((min(max(X, lim(1)), lim(2)) - lim(1))/diff(lim)*(n-1)) + 1; idx(isnan(idx)) = 1;
rgb = reshape(cmap(idx, :), [size(X) 3]); g = repmat(~ok, [1 1 3]); rgb(g) = 0.72;
image(ax, xv, yv, rgb); set(ax, 'YDir', 'normal'); colormap(ax, cmap); clim(ax, lim); colorbar(ax);
end

function [U, D, T] = parts(E, c, P)
Eu = E(:,:,:,1,c,:); Ed = E(:,:,:,2,c,:);
r = @(X) reshape(X, [size(X,1) size(X,2) size(X,3) size(X,6)]);
U = P.*r(sum(abs(Eu).^2, 5)); D = P.*r(sum(abs(Ed).^2, 5)); T = P.*r(sum(abs(Eu + Ed).^2, 5));
end
