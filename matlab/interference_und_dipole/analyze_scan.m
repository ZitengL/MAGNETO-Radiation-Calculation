%% analyze_scan
% Analysis of the undulator + dipole interference scan (run_scan.m). Single electron, 240 MeV, real
% undulator map (centre z = 1 m), dipole (0.1 m, 0.65 T, bend up) at distance D downstream.
% Spectra d2W/(domega dOmega) [J/(rad/s)/sr] per sr of the 11 m sphere centred at the origin.
% Polarisation (observer fields Ex, Ey of the combined run):
%   S0 = |Ex|^2 + |Ey|^2 (+|Ez|^2, negligible), Sx = |Ex|^2, Sy = |Ey|^2,
%   S(+-45) = |Ex +- Ey|^2/2 (linear polariser at +-45 deg), S2 = S(+45) - S(-45) = 2 Re(Ex Ey*),
%   S3 = -2 Im(Ex Ey*) (circular part; sign convention as in Jackson's Stokes parameters for exp(-i w t)).
%   Undulator light is x polarised, dipole light on axis y polarised -> the interference is in S2, S3
%   (and in the 45-degree projections), not in S0. Off axis (0.5 mrad, horizontal) the dipole light has
%   an x component, so S0 shows some interference too.
% Interference term in S0: S0 - S_und - S_dip (undulator alone, dipole alone aligned, see run_scan.m).
% Delay: the fringes are the phase w*tau of Ex Ey*; tau changes with D. Measured change
%   dtau(D) = arg max_t |sum_k z_k exp(i w_k t)|, z_k = phase of (Ex Ey*)(D) * conj((Ex Ey*)(D0)),
%   over the frequencies where |Ex| and |Ey| are both > 10% of their envelope; D0 = 0.5 m.
% Predicted change: arrival-time change at the observer of light leaving the dipole entrance with the
%   electron (t_arr_entrance from run_scan.m, exact geometry at 11 m).
% Outputs: results/interference_scan.mat, results/D0p5_<obs>.txt, figures/*.png/.fig

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
c0 = 299792458; eps0 = 8.85418781e-12; gam = 470.66442480;
U = load(fullfile(res_dir, 'undulator_only.mat'));
lam = U.lambda_nm; w = 2*pi*c0 ./ (lam*1e-9); Nw = numel(lam); obs = U.obs; oname = U.obs_name;
olab = {'on axis (0, 0, 11) m', '\theta_x = 0.5 mrad at 11 m'};
fl = dir(fullfile(res_dir, 'scan', 'D*mm.mat'));
nD = numel(fl); Dv = zeros(1, nD); Ec = complex(zeros(nD, 2, 3, Nw, 'single')); Ed = Ec; tarr = zeros(nD, 2);
for i = 1:nD
    s = load(fullfile(fl(i).folder, fl(i).name)); Dv(i) = s.D; Ec(i,:,:,:) = s.Ew_comb; Ed(i,:,:,:) = s.Ew_dip; tarr(i,:) = s.t_arr_entrance;
end
[Dv, o] = sort(Dv); Ec = Ec(o,:,:,:); Ed = Ed(o,:,:,:); tarr = tarr(o,:);
fprintf('%d D values: %.3f .. %.3f m\n', nD, Dv(1), Dv(end));
xt = [15 20 30 40 50 70 100 150];
nw = round(0.01 / mean(diff(log(lam)))); envf = @(S) movmax(S, 2*nw + 1, 2);
i0 = find(abs(Dv - 0.5) < 1e-9);

A = struct();
for q = 1:2
    pref = eps0*c0*norm(obs{q})^2/pi;
    Ex = double(squeeze(Ec(:,q,1,:))); Ey = double(squeeze(Ec(:,q,2,:))); Ez = double(squeeze(Ec(:,q,3,:)));
    Dx = double(squeeze(Ed(:,q,1,:))); Dy = double(squeeze(Ed(:,q,2,:))); Dz = double(squeeze(Ed(:,q,3,:)));
    Eu = double(squeeze(U.Ew_und(q,:,:)));
    S0 = pref*(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2); Sx = pref*abs(Ex).^2; Sy = pref*abs(Ey).^2;
    Sp45 = pref*abs(Ex + Ey).^2/2; Sm45 = pref*abs(Ex - Ey).^2/2;
    S2 = 2*pref*real(Ex.*conj(Ey)); S3 = -2*pref*imag(Ex.*conj(Ey));
    Su = pref*sum(abs(Eu).^2, 1); Sd = pref*(abs(Dx).^2 + abs(Dy).^2 + abs(Dz).^2);
    Iint = S0 - Su - Sd;                                      % interference term in the total intensity
    cphi = S2 ./ (2*sqrt(Sx.*Sy));                            % cosine of the phase difference between Ex and Ey (fringes)
    Vis  = 2*sqrt(Sx.*Sy) ./ (Sx + Sy);                       % fringe visibility behind a 45-degree polariser
    dec = max(abs(Ex - Eu(1,:) - Dx) ./ envf(abs(Ex)), [], 2); dec_y = max(abs(Ey - Eu(2,:) - Dy) ./ envf(abs(Ey)), [], 2);
    fprintf('\n=== %s\n', oname{q});
    fprintf('decomposition |E_comb - E_und - E_dip|/envelope: max over all D: x %.1e, y %.1e\n', max(dec), max(dec_y));
    fprintf('D = 0.5 m: S_dip/S_und median %.1e; |S0 - S_und - S_dip|/S0 max %.2e (median %.1e); sqrt(S2^2+S3^2)/S0 median %.2f, max %.2f\n', ...
        median(Sd(i0,:)./Su), max(abs(Iint(i0,:))./S0(i0,:)), median(abs(Iint(i0,:))./S0(i0,:)), median(hypot(S2(i0,:), S3(i0,:))./S0(i0,:)), max(hypot(S2(i0,:), S3(i0,:))./S0(i0,:)));
    Ixyz = 2*pref*real([Eu(1,:).*conj(Dx(i0,:)); Eu(2,:).*conj(Dy(i0,:)); Eu(3,:).*conj(Dz(i0,:))]);
    fprintf('           parts of the total-intensity interference / S0, max: x-x %.1e, y-y %.1e, z-z (longitudinal) %.1e\n', max(abs(Ixyz) ./ S0(i0,:), [], 2));

    % ---- delay change relative to D = 0.5 m
    C = Ex.*conj(Ey);
    mask = abs(Ex) > 0.1*envf(abs(Ex)) & abs(Ey) > 0.1*envf(abs(Ey));
    dtau = zeros(1, nD); dtau_pred = tarr(:,q)' - tarr(i0,q);
    for i = 1:nD
        m = mask(i,:) & mask(i0,:); z = C(i,m).*conj(C(i0,m)); z = z./abs(z); wm = w(m);
        tg = dtau_pred(i) + linspace(-2e-15, 2e-15, 2001);   % search +-0.6 um around the prediction
        F = abs(exp(-1i*tg'*wm) * z.');                      % fields ~ exp(-i w t_arrival): a later dipole gives
        dtau(i) = peak_pos(tg, F);                           % Ex Ey* ~ exp(+i w dtau)
    end
    pf = polyfit(Dv, dtau*c0, 1); pp = polyfit(Dv, dtau_pred*c0, 1);
    thd = atan(obs{q}(1) ./ (obs{q}(3) - (1.32004 + Dv + 0.05)));          % angle seen from the dipole
    slope_simple = mean((1/sqrt(1 - 1/gam^2) - cos(thd)));                   % (1/beta - cos theta_d), c dtau/dD
    fprintf('delay change c*dtau vs D: measured slope %.6e, predicted (geometry) %.6e, (1/beta - cos theta_d) %.6e; max |measured - predicted| = %.2e um\n', ...
        pf(1), pp(1), slope_simple, max(abs(dtau - dtau_pred))*c0*1e6);

    % ---- delay spectrum (absolute delay between undulator and dipole light), selected D
    tgrid = linspace(-2e-15, 25e-15, 2701); Dsel = [0.1 0.5 1.0 1.5];
    Dsel = Dsel(ismember(round(Dsel*1e3), round(Dv*1e3))); DS = zeros(numel(Dsel), numel(tgrid)); tpk = zeros(1, numel(Dsel)); tedge = zeros(numel(Dsel), 2);
    dw = abs(gradient(w));
    for j = 1:numel(Dsel)
        i = find(abs(Dv - Dsel(j)) < 1e-9);
        DS(j,:) = abs(exp(-1i*tgrid'*w) * (C(i,:).*dw).');  % t = (arrival of dipole light) - (arrival of undulator light)
        % the undulator light is a wave train (~N lambda_1 long), the dipole light a short flash -> a box;
        % its edges = delays between the flash and the end / start of the train (half maximum), centre = midpoint
        k = find(DS(j,:) > 0.5*max(DS(j,:))); tedge(j,:) = tgrid([k(1) k(end)]); tpk(j) = mean(tedge(j,:));
    end
    fprintf('delay spectrum (box, half maximum) at D = %s m: c*t from %s to %s um; centre %s um (centre differences to D = 0.5: %s um, predicted %s um)\n', ...
        mat2str(Dsel), mat2str(tedge(:,1)'*c0*1e6, 4), mat2str(tedge(:,2)'*c0*1e6, 4), mat2str(tpk*c0*1e6, 4), ...
        mat2str((tpk - tpk(Dsel == 0.5))*c0*1e6, 4), mat2str(interp1(Dv, dtau_pred, Dsel)*c0*1e6, 4));

    A.(oname{q}) = struct('D', Dv, 'S0', S0, 'Sx', Sx, 'Sy', Sy, 'Sp45', Sp45, 'Sm45', Sm45, 'S2', S2, 'S3', S3, ...
        'S_und', Su, 'S_dip', Sd, 'I_int', Iint, 'cos_dphi', cphi, 'visibility45', Vis, 'dtau', dtau, 'dtau_pred', dtau_pred, 'slope_meas', pf(1), 'slope_pred', pp(1), ...
        'delay_t', tgrid, 'delay_spec', DS, 'delay_D', Dsel, 'delay_peak', tpk);

    % ---- txt for D = 0.5 m
    fid = fopen(fullfile(res_dir, ['D0p5_' oname{q} '.txt']), 'wt');
    fprintf(fid, '%% Single electron 240 MeV, real undulator (centre z = 1 m) + dipole (0.1 m, 0.65 T, bend up) at D = 0.5 m downstream; observer %s\n', strrep(olab{q}, '\theta_x', 'theta_x'));
    fprintf(fid, '%% d2W/domega dOmega [J/(rad/s)/sr] per sr of the 11 m sphere centred at the origin\n');
    fprintf(fid, '%% Columns: 1 lambda [nm]  2 omega [rad/s]  3 S0 total  4 Sx  5 Sy  6 S(+45 deg)  7 S(-45 deg)  8 S2  9 S3  10 undulator alone  11 dipole alone  12 S0 - und - dip\n');
    fprintf(fid, ['%-14.8f %-14.6e' repmat(' %-14.6e', 1, 10) '\n'], [lam; w; S0(i0,:); Sx(i0,:); Sy(i0,:); Sp45(i0,:); Sm45(i0,:); S2(i0,:); S3(i0,:); Su; Sd(i0,:); Iint(i0,:)]);
    fclose(fid);

    % ---- figure 1: D = 0.5 m spectra
    f = figure('Position', [100 100 1000 900]);
    subplot(3,1,1);
    loglog(lam, S0(i0,:), 'k-', lam, Su, 'b--', lam, Sd(i0,:), 'r-', lam, Su + Sd(i0,:), 'g:', 'LineWidth', 1.1); grid on; xlim([15 150]); set(gca, 'XTick', xt);
    legend('undulator + dipole (total)', 'undulator alone', 'dipole alone', 'incoherent sum', 'Location', 'southeast');
    ylabel('d^2W/d\omega d\Omega'); title(['Undulator + dipole at D = 0.5 m, single electron 240 MeV, ' olab{q}]);
    subplot(3,1,2);
    loglog(lam, Sx(i0,:), 'b-', lam, Sy(i0,:), 'r-', lam, Sp45(i0,:), 'k-', lam, Sm45(i0,:), 'm-'); grid on; xlim([15 150]); set(gca, 'XTick', xt);
    legend('x polarisation', 'y polarisation', '+45 deg polariser', '-45 deg polariser', 'Location', 'southeast'); ylabel('d^2W/d\omega d\Omega');
    subplot(3,1,3);
    semilogx(lam, cphi(i0,:), 'k-', lam, Vis(i0,:), 'r-', lam, Iint(i0,:)./S0(i0,:), 'b-'); grid on; xlim([15 150]); set(gca, 'XTick', xt); ylim([-1.05 1.05]);
    legend('cos\Delta\phi = S_2 / (2(S_x S_y)^{1/2})', 'visibility with a 45^o polariser, 2(S_x S_y)^{1/2}/(S_x+S_y)', '(S_0 - S_{und} - S_{dip}) / S_0', 'Location', 'southoutside', 'NumColumns', 3);
    xlabel('\lambda [nm]'); ylabel('normalised');
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(f, fullfile(fig_dir, ['D0p5_spectra_' oname{q} '.png']), 'Resolution', 200); savefig(f, fullfile(fig_dir, ['D0p5_spectra_' oname{q} '.fig'])); close(f);

    % ---- figure 2: (lambda, D) maps, fine scan
    fine = Dv >= 0.5 - 1e-9 & Dv <= 0.7 + 1e-9;
    f = figure('Position', [100 100 1000 800]);
    subplot(2,1,1);
    imagesc(log10(lam), Dv(fine), cphi(fine,:)); axis xy; colormap(gca, redblue()); caxis([-1 1]); cb = colorbar; cb.Label.String = 'cos\Delta\phi = S_2/(2(S_xS_y)^{1/2})';
    set(gca, 'XTick', log10(xt), 'XTickLabel', xt); ylabel('D [m]'); title(['fringes: phase difference between E_x (undulator) and E_y (dipole), ' olab{q}]);
    subplot(2,1,2);
    imagesc(log10(lam), Dv(fine), Iint(fine,:)./S0(fine,:)); axis xy; colormap(gca, redblue()); cb = colorbar; cb.Label.String = '(S_0 - S_{und} - S_{dip})/S_0';
    v = max(abs(Iint(fine,:)./S0(fine,:)), [], 'all'); caxis([-v v]);
    set(gca, 'XTick', log10(xt), 'XTickLabel', xt); xlabel('\lambda [nm]'); ylabel('D [m]'); title('interference term in the total intensity');
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(f, fullfile(fig_dir, ['map_lambda_D_' oname{q} '.png']), 'Resolution', 200); savefig(f, fullfile(fig_dir, ['map_lambda_D_' oname{q} '.fig'])); close(f);
    % zooms: the full-band map is too dense to render at short wavelengths (moire); narrow windows
    win = [15 16.5; 40 45; 100 150];
    f = figure('Position', [100 100 1200 450]);
    for j = 1:3
        subplot(1,3,j); sel = lam >= win(j,1) & lam <= win(j,2);
        imagesc(lam(sel), Dv(fine), cphi(fine, sel)); axis xy; colormap(gca, redblue()); caxis([-1 1]);
        xlabel('\lambda [nm]'); if j == 1; ylabel('D [m]'); end
        title(sprintf('%g-%g nm', win(j,:)));
        fr = 2*gam^2*win(j,1)*1e-9 / (1 + (gam*mean(thd))^2);
        text(0.02, 0.02, sprintf('1 fringe per %.1f mm in D at %g nm', fr*1e3, win(j,1)), 'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'w', 'VerticalAlignment', 'bottom');
    end
    cb = colorbar; cb.Label.String = 'cos\Delta\phi';
    sgtitle(['fringes (cos\Delta\phi), zoomed, 1 mm steps in D, ' olab{q}]);
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(f, fullfile(fig_dir, ['map_zoom_' oname{q} '.png']), 'Resolution', 200); savefig(f, fullfile(fig_dir, ['map_zoom_' oname{q} '.fig'])); close(f);

    % ---- figure 3: delay change vs D, and delay spectrum
    f = figure('Position', [100 100 1000 800]);
    subplot(2,1,1);
    plot(Dv, dtau*c0*1e6, 'ko', Dv, dtau_pred*c0*1e6, 'r-', 'MarkerSize', 3); grid on;
    legend('measured from the fringes (phase of E_x E_y^*)', 'predicted: arrival time of light from the dipole entrance', 'Location', 'northwest');
    xlabel('D [m]'); ylabel('c \Delta\tau relative to D = 0.5 m [\mum]');
    title(sprintf('%s: slope measured %.4e, predicted %.4e (\\approx (1+\\gamma^2\\theta_d^2)/(2\\gamma^2))', olab{q}, pf(1), pp(1)));
    subplot(2,1,2);
    plot(tgrid*c0*1e6, DS ./ max(DS, [], 2), 'LineWidth', 1.1); grid on; xlim([tgrid(1) tgrid(end)]*c0*1e6);
    legend(arrayfun(@(d) sprintf('D = %.1f m', d), Dsel, 'UniformOutput', false), 'Location', 'northeast');
    xlabel('c t [\mum]'); ylabel('|\Sigma E_x E_y^* e^{i\omega t} \Delta\omega| (normalised)'); title('delay spectrum: delay between undulator and dipole light');
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(f, fullfile(fig_dir, ['delay_vs_D_' oname{q} '.png']), 'Resolution', 200); savefig(f, fullfile(fig_dir, ['delay_vs_D_' oname{q} '.fig'])); close(f);
end
A.lambda_nm = lam; A.omega = w; A.note = 'per sr of the 11 m sphere; S2 = 2 Re(Ex Ey*), S3 = -2 Im(Ex Ey*); dtau relative to D = 0.5 m';
save(fullfile(res_dir, 'interference_scan.mat'), 'A', '-v7.3');

function tp = peak_pos(t, F)
% position of the maximum of F(t) on a uniform grid t, refined with a parabola in index units
[~, j] = max(F); j = min(max(j, 2), numel(t)-1);
d = 0.5*(F(j-1) - F(j+1)) / (F(j-1) - 2*F(j) + F(j+1));
tp = t(j) + d*(t(2) - t(1));
end

function cm = redblue()
x = linspace(0, 1, 128)';
cm = [[x; ones(128,1)], [x; flipud(x)], [ones(128,1); flipud(x)]];
end
