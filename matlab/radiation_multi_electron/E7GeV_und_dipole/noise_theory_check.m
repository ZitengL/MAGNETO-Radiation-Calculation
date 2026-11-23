%% noise_theory_check
% Does the noise in the 5e4 run behave as predicted for macro-particle (finite sampling) noise?  y = -2.0 mm, x = 0, z = 5.9 m.
% Theory (N_m macro-particles with independent random positions; a_j = field of particle j, phi_j its phase):
%   <|sum_j a_j e^{i phi_j}|^2> = sum_j |a_j|^2  +  N_m (N_m - 1) |<a e^{i phi}>|^2  =  self terms  +  coherent part C,
%   and for large N_m the noise term is a complex Gaussian => for fixed C the spectrum at one frequency is
%   |sqrt(C) + n|^2, n ~ CN(0, self) (non-central chi-square with 2 degrees of freedom per field component).
%   With the time shift, C(sigma_z) = C(10 nm) exp(-k^2 (sigma_z^2 - (10 nm)^2)) (Gaussian X3 independent of the rest).
% Tests:
%   1) spectrum: S_sim(sigma_z) vs the prediction C_pred + self with its 16-84 % band (Monte Carlo of |sqrt(C) + n|^2);
%      C(10 nm) = S_sim(10 nm) - self (at 10 nm S >> self over the whole band).
%   2) noise energy: from the 25 chunks, A = sum_c |E_c|^2 = self + C/25 and B = |sum_c E_c|^2 = self + C, so the
%      noise level can be measured independently of the self terms: self_est = (25 A - B)/24; compare W(self_est) with W(self).
%   3) distribution: in the deep-noise region (chunk coherent part < 1 % of the chunk self terms) |E_c|^2 / self_c
%      should be exponential with mean 1 and standard deviation 1.
% x polarisation (what the polariser passes). Output: figures/noise_theory_y-2.0mm.png/.fig, results/noise_theory_y-2.0mm.txt

clear; close all;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; eps0 = 8.85418781e-12; yo = -2.0;
R = load(fullfile(res_dir, 'combined.mat')); R = R.R;
om = R.omega; lnm = R.lambda*1e9; k = om/c0; sz = R.sigma_z; ns = numel(sz); nch = size(R.Ec, 1); nw = numel(om);
o = find(abs(1e3*R.obs(:,2) - yo) < 1e-9); P = eps0*c0*norm(R.obs(o,:))^2/pi;
Wn = @(v) -trapz(om, v);
Ec = squeeze(R.Ec(:,:,o,1,1,:) + R.Ec(:,:,o,2,1,:));               % chunk x sz x omega, x component, und + dip
selfc = P*squeeze(R.Ic(:,o,3,1,:));                                 % chunk x omega, self terms of each chunk
self = sum(selfc, 1);                                               % 1 x omega
B = P*squeeze(abs(sum(Ec, 1)).^2);                                  % sz x omega (full sum = S_sim)
A = P*squeeze(sum(abs(Ec).^2, 1));                                  % sz x omega
C10 = B(1,:) - self;
Cp = C10.*exp(-k.^2.*(sz(:).^2 - sz(1)^2));                         % sz x omega
rng(1); nmc = 4000; n = (randn(nmc, nw) + 1i*randn(nmc, nw))/sqrt(2);
lo = zeros(ns, nw); hi = lo;
for s = 1:ns
    X = abs(sqrt(Cp(s,:)) + sqrt(self).*n).^2;
    X = sort(X, 1); lo(s,:) = X(round(0.16*nmc),:); hi(s,:) = X(round(0.84*nmc),:);
end
self_est = (nch*A - B)/(nch - 1);

% 3) distribution in the deep-noise region
r = [];
for s = 1:ns
    m = Cp(s,:)/nch^2 < 0.01*self/nch;                              % chunk coherent part < 1 % of chunk self terms
    r = [r; reshape(P*abs(squeeze(Ec(:,s,m))).^2./selfc(:,m), [], 1)]; %#ok<AGROW>
end

fid = fopen(fullfile(res_dir, sprintf('noise_theory_y%.1fmm.txt', yo)), 'w');
fprintf(fid, '%% x = 0, y = %.1f mm, z = 5.9 m, x polarisation, 5e4 macro-particles in 25 chunks; energies in 50-300 nm [J/sr]\n', yo);
fprintf(fid, '%% predicted noise energy W(self) = %.4e J/sr (independent of sigma_z)\n', Wn(self));
fprintf(fid, '%% sigma_z_nm  W_sim  W_pred(=C_pred+self)  W_C_pred  W_self_est(from chunks)  W_self_est/W_self  lambda_cross_nm(C_pred=self)\n');
for s = 1:ns
    d = Cp(s,:) - self; ic = find(d > 0, 1);                        % lambda ascending: first lambda where C_pred > self
    if isempty(ic); lc = nan; elseif ic == 1; lc = lnm(1); else; lc = lnm(ic); end
    fprintf(fid, '%3.0f  %.4e  %.4e  %.4e  %.4e  %.3f  %6.1f\n', 1e9*sz(s), Wn(B(s,:)), Wn(Cp(s,:) + self), Wn(Cp(s,:)), Wn(self_est(s,:)), Wn(self_est(s,:))/Wn(self), lc);
end
fprintf(fid, '%% deep-noise region, |E_c|^2/self_c over chunks, sigma_z and wavelengths: %d values, mean %.3f, s.d. %.3f, median %.3f (exponential: 1, 1, 0.693)\n', numel(r), mean(r), std(r), median(r));
fprintf(fid, '%% fraction of S_sim inside the predicted 16-84 %% band, per sigma_z (all wavelengths):');
for s = 1:ns; fprintf(fid, ' %.2f', mean(B(s,:) >= lo(s,:) & B(s,:) <= hi(s,:))); end
fprintf(fid, '\n');
fclose(fid); type(fullfile(res_dir, sprintf('noise_theory_y%.1fmm.txt', yo)));

f = figure('Position', [30 30 1500 1150], 'Visible', 'off');
tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, {sprintf('x = 0, y = %.1f mm, z = 5.9 m, x polariser: simulation vs macro-particle noise theory', yo), ...
    'prediction = C_{pred} + self terms, band = 16-84 % of |\surdC_{pred} + n|^2, n ~ CN(0, self);  C_{pred} = C(10 nm) exp(-k^2(\sigma_z^2 - (10 nm)^2))'});
tk = [50 70 100 150 200 300];
for s = 1:ns
    ax(s) = nexttile; hold on; %#ok<SAGROW>
    fill([lnm fliplr(lnm)], [max(lo(s,:), 1e-30) fliplr(hi(s,:))], 'r', 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'DisplayName', 'theory 16-84 %');
    plot(lnm, Cp(s,:) + self, 'r--', 'LineWidth', 1.2, 'DisplayName', 'theory mean C_{pred} + self');
    plot(lnm, Cp(s,:), 'b:', 'LineWidth', 1.2, 'DisplayName', 'C_{pred} (coherent only)');
    plot(lnm, self, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4, 'DisplayName', 'self terms');
    plot(lnm, self_est(s,:), '-', 'Color', [0 0.6 0], 'LineWidth', 0.8, 'DisplayName', 'noise measured from chunks');
    plot(lnm, B(s,:), 'k-', 'LineWidth', 1.6, 'DisplayName', 'simulation S_{tot,x}');
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', tk, 'Layer', 'top'); grid on; box on; xlim([50 300]);
    xlabel('\lambda [nm]'); ylabel('d^2W/d\omegad\Omega [J/(rad/s)/sr]'); title(sprintf('\\sigma_z = %g nm', 1e9*sz(s)));
    if s == 1; legend('Location', 'southeast', 'FontSize', 7); end
end
set(ax, 'YLim', [1e-20 1e-14]);
saveas(f, fullfile(fig_dir, sprintf('noise_theory_y%.1fmm.png', yo))); savefig(f, fullfile(fig_dir, sprintf('noise_theory_y%.1fmm.fig', yo)));
