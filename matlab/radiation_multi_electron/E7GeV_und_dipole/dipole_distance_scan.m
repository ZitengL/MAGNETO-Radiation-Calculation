%% dipole_distance_scan
% User's question (2026-11-05): does moving the dipole change the wavelength of the undulator-dipole interference
% notch, can the two add up, and what does it do to the bunch-length measurement?
% Single ideal electron (undulator turned over, result_sm_dz60.gdf), observer x = 0, y = -2.0 mm, z = 5.9 m, no polariser.
% Moving the dipole downstream by dL: between the magnets the electron drifts in a field-free straight line, and the
% rectmagnet field depends only on z (uniform across the beam), so the dipole part of the trajectory is the old one
% shifted by v0*dt in space and dt = dL/v0z in time (v0 = drift velocity before the dipole). This is exact for the ideal
% electron (the small Bz ~ y dB/dz fringe term is neglected for the transverse offset x' dL ~ 0.5 um). Only dL >= 0:
% moving the dipole upstream would put its fringe field into the undulator (user's constraint).
% Bunch-length part: a MODEL, not the 5e4 simulation: coherent spectrum N^2 S1(lambda) exp(-k^2 sz^2) exp(-k^2 s^2 th^2),
% s = 20 um, th = angle of the observer from the undulator centre (as in the positions-only model of analyze_combined.m).
% Output: figures/dipole_distance_scan.png/.fig, results/dipole_distance_scan.txt

clear; close all;
this_dir = fileparts(mfilename('fullpath')); root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
c0 = 299792458; eps0 = 8.85418781e-12; N = 100e-12/1.602176634e-19;
lam = logspace(log10(50), log10(300), 300)*1e-9; om = 2*pi*c0./lam; k = om/c0; lnm = lam*1e9;
obs = [0 -2e-3 5.9]; pref = eps0*c0*norm(obs)^2/pi; Wn = @(v) -trapz(om, v);
tr = read_gpt_trajectory(fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole', 'ideal_electron', 'result_sm_dz60.gdf'));
su = sub(tr, tr.r(3,:) < 0.63); sd0 = sub(tr, tr.r(3,:) >= 0.63);
v0 = c0*sd0.b(:,1);                                                  % drift velocity before the dipole
Eu = emission_time_spectrum(su, obs, lam, 'fB', norm(obs)/c0).Ew;   % 3 x nw
Su = pref*sum(abs(Eu).^2, 1);
dLs = [0 0.05 0.1 0.2 0.3 0.5 0.75 1.0 1.5 2.0];
szs = [10 20 30 40 50]*1e-9; th = norm(obs(1:2))/(obs(3) - 0.4); Ft = exp(-k.^2*(20e-6)^2*th^2);
fid = fopen(fullfile(this_dir, 'results', 'dipole_distance_scan.txt'), 'w');
fprintf(fid, '%% single ideal electron, x = 0, y = -2.0 mm, z = 5.9 m, no polariser; dipole moved downstream by dL (dipole centre 1.3144 m + dL)\n');
fprintf(fid, '%% W per electron in 50-300 nm [J/sr]; notch = wavelengths of local minima of S_tot below 0.3 (S_und + S_dip)\n');
fprintf(fid, '%% model W_coh(sz) = int N^2 S1 exp(-k^2 sz^2) Ft (positions only), sz = 10 20 30 40 50 nm, and ratio 10/50\n');
fprintf(fid, '%% dL_m W_und W_dip W_tot W_tot/(W_und+W_dip) y_beam_at_5.9m_mm | Wcoh_tot(10..50 nm) | ratio_tot ratio_und | notches_nm\n');
St = zeros(numel(dLs), numel(lam)); Sd = St;
for i = 1:numel(dLs)
    dt = dLs(i)/v0(3); s = sd0; s.t = s.t + dt; s.r = s.r + v0*dt;
    Ed = emission_time_spectrum(s, obs, lam, 'fB', norm(obs)/c0).Ew;
    Sd(i,:) = pref*sum(abs(Ed).^2, 1); St(i,:) = pref*sum(abs(Eu + Ed).^2, 1);
    yb = s.r(2,end) + s.b(2,end)/s.b(3,end)*(5.9 - s.r(3,end));
    q = St(i,:)./(Su + Sd(i,:)); mn = find(islocalmin(St(i,:)) & q < 0.3);
    Wc = arrayfun(@(z) Wn(N^2*St(i,:).*exp(-k.^2*z^2).*Ft), szs); Wcu = arrayfun(@(z) Wn(N^2*Su.*exp(-k.^2*z^2).*Ft), szs);
    fprintf(fid, '%4.2f %.4e %.4e %.4e %.3f %6.2f | %s | %6.1f %6.1f | %s\n', dLs(i), Wn(Su), Wn(Sd(i,:)), Wn(St(i,:)), Wn(St(i,:))/(Wn(Su) + Wn(Sd(i,:))), ...
        1e3*yb, sprintf('%.3e ', Wc), Wc(1)/Wc(end), Wcu(1)/Wcu(end), sprintf('%.0f ', lnm(mn)));
end
fclose(fid); type(fullfile(this_dir, 'results', 'dipole_distance_scan.txt'));

f = figure('Position', [30 30 1500 600], 'Visible', 'off');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'single ideal electron, x = 0, y = -2.0 mm, z = 5.9 m, no polariser: dipole moved downstream by \DeltaL');
nexttile; hold on; sel = [1 3 5 6 8 10]; cl = turbo(numel(sel) + 1);
plot(lnm, Su, 'b-', 'LineWidth', 2, 'DisplayName', 'undulator part');
plot(lnm, Sd(1,:), 'r-', 'LineWidth', 2, 'DisplayName', 'dipole part (\DeltaL = 0)');
for j = 1:numel(sel)
    plot(lnm, St(sel(j),:), '-', 'Color', cl(j,:), 'LineWidth', 1.2, 'DisplayName', sprintf('total, \\DeltaL = %.2f m', dLs(sel(j))));
end
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', [50 70 100 150 200 300]); grid on; box on; xlim([50 300]);
xlabel('\lambda [nm]'); ylabel('d^2W/d\omegad\Omega per electron [J/(rad/s)/sr]'); legend('Location', 'southwest', 'FontSize', 8);
title('spectra');
nexttile; hold on;
for j = 1:numel(sel)
    plot(lnm, St(sel(j),:)./(Su + Sd(sel(j),:)), '-', 'Color', cl(j,:), 'LineWidth', 1.2, 'DisplayName', sprintf('\\DeltaL = %.2f m', dLs(sel(j))));
end
yline(1, 'k:', 'HandleVisibility', 'off'); set(gca, 'XScale', 'log', 'XTick', [50 70 100 150 200 300]); grid on; box on; xlim([50 300]); ylim([0 2]);
xlabel('\lambda [nm]'); ylabel('S_{tot} / (S_{und} + S_{dip})'); title('interference (1 = none, 2 = fully constructive, 0 = fully destructive)');
legend('Location', 'southwest', 'FontSize', 8);
saveas(f, fullfile(this_dir, 'figures', 'dipole_distance_scan.png')); savefig(f, fullfile(this_dir, 'figures', 'dipole_distance_scan.fig'));

function t = sub(t, m)
k = find(m); kf = k(any(t.fB(:, k) ~= 0, 1)); k = max(k(1), kf(1)-1) : min(k(end), kf(end)+1);
for f = {'t', 'G'};             t.(f{1}) = t.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; t.(f{1}) = t.(f{1})(:, k); end
end
