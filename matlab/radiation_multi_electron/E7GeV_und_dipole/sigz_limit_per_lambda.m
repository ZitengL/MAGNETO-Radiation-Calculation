%% sigz_limit_per_lambda
% User question (2026-11-06): y = -2.0 mm does not cover sigma_z = 10-1000 nm; is there another point that does?
% (a) Bound valid for ANY observation point and polariser (z Gaussian, independent of the other coordinates, no SC):
%     coherent / incoherent = N |F(k)|^2 * eta,  eta = |sum_j <E_j>|^2 / (N sum_j <|E_j|^2>) <= 1 (Cauchy-Schwarz),
%     |F|^2 = exp(-k^2 sz^2), N = 100 pC / e. So coherent = r x incoherent at sz <= sqrt(ln(N/r))/k.
% (b) At y = -2.0 mm (x pol., data of run_large_sigz.m): per wavelength
%     R0(lambda) = (|sum E0|^2 - sum |E0|^2) nmacro / sum |E0|^2 = N eta,  sz_r(lambda) = sqrt(ln(R0/r))/k,
%     i.e. a spectrally resolved detector at that wavelength. Jackknife over the 25 chunks.
% Output: results/sigz_limit_per_lambda.txt, figures/sigz_limit_per_lambda.png/.fig
clear; close all;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; N = 100e-12/1.602176634e-19;
R = load(fullfile(res_dir, 'large_sigz.mat')); R = R.R; nch = numel(R.done);
lam = R.lambda; k = R.omega/c0;
E0c = squeeze(R.E0(:,1,1,:) + R.E0(:,2,1,:)); Ic = squeeze(R.I(:,3,1,:));
R0f = @(E0, I, nm) (abs(E0).^2 - I).*nm./I;
rr = [10 1 0.1];
lim = @(R0, r) real(sqrt(max(log(R0/r), 0)))./k;
R0 = R0f(sum(E0c, 1), sum(Ic, 1), R.nmacro); eta = R0/N;
f = nch/(nch-1); szj = zeros(nch, numel(rr), numel(lam));
for c = 1:nch
    m = setdiff(1:nch, c); R0j = R0f(f*sum(E0c(m,:), 1), f^2*sum(Ic(m,:), 1), R.nmacro*f);
    for i = 1:numel(rr); szj(c,i,:) = lim(R0j, rr(i)); end
end
sd = squeeze(sqrt((nch-1)/nch*sum((szj - mean(szj, 1)).^2, 1)));
fid = fopen(fullfile(res_dir, 'sigz_limit_per_lambda.txt'), 'w');
fprintf(fid, '%% sigma_z [nm] where coherent = r x incoherent at a single wavelength. N = %.3e\n', N);
fprintf(fid, '%% bound (any point, eta = 1) | y = -2.0 mm x pol. (+- jackknife), eta = coherent efficiency\n');
fprintf(fid, '%% lambda_nm | bound r=10 1 0.1 | y-2mm r=10 sd r=1 sd r=0.1 sd | eta\n');
for l = [100 200 280 300 400 500 600 700 800 900 1000]
    [~, i] = min(abs(lam - l*1e-9));
    fprintf(fid, '%6.0f | %5.0f %5.0f %5.0f |', lam(i)*1e9, 1e9*sqrt(log(N./rr))/k(i));
    for j = 1:3; L = lim(R0, rr(j)); fprintf(fid, ' %5.0f %3.0f', 1e9*L(i), 1e9*sd(j,i)); end
    fprintf(fid, ' | %.2e\n', eta(i));
end
for s = [1000e-9]
    fprintf(fid, '%% wavelength needed for sigma_z = %.0f nm: eta = 1: coherent = 1 x inc. at %.2f um, 0.1 x at %.2f um\n', s*1e9, ...
        2*pi*s/sqrt(log(N))*1e6, 2*pi*s/sqrt(log(N/0.1))*1e6);
    fprintf(fid, '%%   if eta stayed at its 1 um value at y = -2 mm (%.1e; an assumption, eta decreases with lambda here): coherent = 1 x inc. at %.2f um, 0.1 x at %.2f um\n', ...
        eta(end), 2*pi*s/sqrt(log(N*eta(end)))*1e6, 2*pi*s/sqrt(log(N*eta(end)/0.1))*1e6);
end
fclose(fid); type(fullfile(res_dir, 'sigz_limit_per_lambda.txt'));
fh = figure('Position', [100 100 1000 420], 'Visible', 'off'); col = lines(3);
subplot(1,2,1);
for j = 1:3
    plot(lam*1e9, 1e9*sqrt(log(N/rr(j)))./k, ':', 'Color', col(j,:), 'LineWidth', 1.5, 'DisplayName', sprintf('bound, r = %g', rr(j))); hold on;
    plot(lam*1e9, 1e9*lim(R0, rr(j)), '-', 'Color', col(j,:), 'DisplayName', sprintf('y = -2 mm, r = %g', rr(j)));
end
yline(1000, 'k--', 'HandleVisibility', 'off'); xline(280, 'k:', 'HandleVisibility', 'off'); xline(300, 'k:', 'HandleVisibility', 'off');
grid on; xlabel('\lambda [nm]'); ylabel('\sigma_z where coherent = r \times incoherent [nm]'); legend('Location', 'northwest', 'FontSize', 8);
title('single wavelength'); ylim([0 1100]);
subplot(1,2,2); semilogy(lam*1e9, eta, 'k'); grid on; xlabel('\lambda [nm]'); ylabel('\eta = coherent efficiency');
title('y = -2.0 mm, x pol. (\eta \leq 1 for any point)');
sgtitle('Largest measurable \sigma_z per wavelength (100 pC, Gaussian z)');
savefig(fh, fullfile(fig_dir, 'sigz_limit_per_lambda.fig')); exportgraphics(fh, fullfile(fig_dir, 'sigz_limit_per_lambda.png'), 'Resolution', 150);
