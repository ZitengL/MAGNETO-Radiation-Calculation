%% analyze_large_sigz
% Results of run_large_sigz.m: x = 0, y = -2.0 mm, z = 5.9 m, x polariser (x component, undulator + dipole),
% sigma_z = 10:10:1000 nm, SLAC's two detectors 50-280 nm (D1) and 300-1000 nm (D2).
% EXPECTED-VALUE METHOD (approved by the user 2026-11-05): the initial z_j are independent of the other coordinates
% and Gaussian, and without space charge E_j(sz) = E0_j exp(i w sz X3_j / v_zj). Averaging over the z_j:
%   <|sum_j E_j|^2> = |F|^2 (|sum_j E0_j|^2 - sum_j |E0_j|^2) + sum_j |E0_j|^2,   |F|^2 = exp(-(w sz / c)^2)
% For the real beam (N = 6.2e8 electrons, nmacro electrons per macro-particle) the incoherent part is
% sum_j |E_j|^2 / nmacro, so  S = pref [ |F|^2 (|sum E0|^2 - sum |E0|^2) + sum |E0|^2 / nmacro ].
% This removes the longitudinal macro-particle noise; the transverse sampling error is the jackknife over the 25 chunks.
% Checks: (1) one macro-particle realisation at 10-100 nm (Er) vs the formula with the macro-particle self term;
%         (2) Er at 10-50 nm vs the earlier run (combined.mat, 50-300 nm); (3) lambda grid: every 2nd point.
% Output: results/large_sigz_summary.txt/.mat, results/large_sigz_spectra.txt, figures/large_sigz_*.png/.fig

clear; close all;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; eps0 = 8.85418781e-12;
R = load(fullfile(res_dir, 'large_sigz.mat')); R = R.R;
dn = find(R.done); nch = numel(dn); if nch < numel(R.done); fprintf('WARNING: only %d chunks\n', nch); end
lam = R.lambda; om = R.omega; k = om/c0; pref = eps0*c0*norm(R.obs)^2/pi; nm = R.nmacro*numel(R.done)/nch;
sc = numel(R.done)/nch;                                            % keep 100 pC if chunks are missing
E0c = squeeze(R.E0(dn,1,1,:) + R.E0(dn,2,1,:))*sc;                 % chunk x omega, x component, und + dip
Ic  = squeeze(R.I(dn,3,1,:))*sc^2;
sz = (10:10:1000)*1e-9; ns = numel(sz);
bands = [50e-9 280e-9; 300e-9 1000e-9; 50e-9 300e-9]; bn = {'D1 50-280 nm', 'D2 300-1000 nm', '50-300 nm'}; nb = 3;
Wb = @(S, b) -trapz(om(lam >= bands(b,1) & lam <= bands(b,2)), S(lam >= bands(b,1) & lam <= bands(b,2)));
spec = @(E0, I, nmac, s) pref*(exp(-(k*s).^2).*(abs(E0).^2 - I) + I/nmac);

E0 = sum(E0c, 1); I = sum(Ic, 1);
W = zeros(nb, ns); Wcoh = zeros(nb, ns); Winc = zeros(nb, 1); Sall = zeros(ns, numel(lam));
for b = 1:nb; Winc(b) = Wb(pref*I/nm, b); end
for s = 1:ns
    Sall(s,:) = spec(E0, I, nm, sz(s));
    for b = 1:nb; W(b,s) = Wb(Sall(s,:), b); Wcoh(b,s) = W(b,s) - Winc(b); end
end
sens = zeros(nb, ns); for b = 1:nb; sens(b,:) = -gradient(log(W(b,:)), sz)*1e-9; end   % per nm

% jackknife (drop a chunk, rescale to 100 pC: fields x f, self terms x f^2, nmacro x f)
f = nch/(nch-1); Wj = zeros(nch, nb, ns); sj = zeros(nch, nb, ns);
for c = 1:nch
    m = setdiff(1:nch, c); E0j = f*sum(E0c(m,:), 1); Ij = f^2*sum(Ic(m,:), 1);
    for s = 1:ns
        Sj = spec(E0j, Ij, nm*f, sz(s));
        for b = 1:nb; Wj(c,b,s) = Wb(Sj, b); end
    end
    for b = 1:nb; sj(c,b,:) = -gradient(log(squeeze(Wj(c,b,:))'), sz)*1e-9; end
end
jsd = @(A) squeeze(sqrt((nch-1)/nch*sum((A - mean(A, 1)).^2, 1)));
Wsd = jsd(Wj); sens_sd = jsd(sj);

% lambda grid check: every 2nd point
lam2 = lam(1:2:end); om2 = om(1:2:end);
gridchk = zeros(nb, ns);
for b = 1:nb
    m2 = lam2 >= bands(b,1) & lam2 <= bands(b,2);
    for s = 1:ns
        S2 = Sall(s, 1:2:end); gridchk(b,s) = -trapz(om2(m2), S2(m2))/W(b,s) - 1;
    end
end

% check 1: one realisation (Er) vs expected value with the macro-particle self term, 10-100 nm
fid = fopen(fullfile(res_dir, 'large_sigz_summary.txt'), 'w');
fprintf(fid, '%% x = 0, y = -2.0 mm, z = 5.9 m, x polarisation (undulator + dipole), 100 pC, expected-value method\n');
fprintf(fid, '%% CHECK 1: realisation with 5e4 macro-particles vs expected value incl. macro-particle self term (macro self = %.3e J/sr in 50-300 nm)\n', Wb(pref*I, 3));
for s = 1:numel(R.sz_real)
    Er = sum(squeeze(R.Er(dn,s,1,1,:) + R.Er(dn,s,2,1,:)), 1)*sc;
    for b = [1 2]
        a = Wb(pref*abs(Er).^2, b); e = Wb(spec(E0, I, 1, R.sz_real(s)), b);
        fprintf(fid, '%%   sz %3.0f nm, %s: realisation %.4e, expected %.4e, diff/self %+.2f\n', R.sz_real(s)*1e9, bn{b}, a, e, (a - e)/Wb(pref*I, b));
    end
end
% check 2: realisation vs the earlier run combined.mat (same particles and GPT, y = -2 mm is observer 3), 50-300 nm
C = load(fullfile(res_dir, 'combined.mat')); C = C.R; o = find(abs(C.obs(:,2) + 2e-3) < 1e-9);
for s = 1:5                                                        % this run 10:10:50 nm = combined.mat 10:5:50 nm (2s-1)
    cs_ = 2*s - 1;
    Eo = sum(squeeze(C.Ec(dn,cs_,o,1,1,:) + C.Ec(dn,cs_,o,2,1,:)), 1);
    Wo = -trapz(C.omega, pref*abs(Eo).^2);
    Er = sum(squeeze(R.Er(dn,s,1,1,:) + R.Er(dn,s,2,1,:)), 1);
    fprintf(fid, '%% CHECK 2 (same chunks): sz %2.0f nm, 50-300 nm: this run %.5e, combined.mat %.5e, rel %+.2e\n', C.sigma_z(cs_)*1e9, Wb(pref*abs(Er).^2, 3), Wo, Wb(pref*abs(Er).^2, 3)/Wo - 1);
end
fprintf(fid, '%% CHECK 3: lambda grid (every 2nd point) max |rel. change| of W: D1 %.1e, D2 %.1e\n', max(abs(gridchk(1,:))), max(abs(gridchk(2,:))));
fprintf(fid, '%% incoherent level (real beam): D1 %.3e, D2 %.3e J/sr\n', Winc(1), Winc(2));
for b = [1 2]
    for r = [10 1 0.1]
        i = find(Wcoh(b,:) < r*Winc(b), 1);
        if isempty(i); fprintf(fid, '%% %s: coherent > %g x incoherent up to 1000 nm\n', bn{b}, r);
        else; fprintf(fid, '%% %s: coherent < %g x incoherent from sigma_z = %.0f nm\n', bn{b}, r, sz(i)*1e9); end
    end
end
fprintf(fid, '%% sigma_z_nm | W_D1 sd coh/inc sens_%%/nm sens_sd | W_D2 sd coh/inc sens_%%/nm sens_sd  (W in J/sr)\n');
for s = 1:ns
    fprintf(fid, '%5.0f', sz(s)*1e9);
    for b = [1 2]; fprintf(fid, ' | %.4e %.2e %.3e %7.3f %.3f', W(b,s), Wsd(b,s), Wcoh(b,s)/Winc(b), 100*sens(b,s), 100*sens_sd(b,s)); end
    fprintf(fid, '\n');
end
fclose(fid); type(fullfile(res_dir, 'large_sigz_summary.txt'));
sel = [10 20 30 50 70 100 150 200 300 400 500 600 700 800 1000];
[~, is] = ismember(round(sel), round(sz*1e9));
T = [lam(:)*1e9, Sall(is,:)'];
fid = fopen(fullfile(res_dir, 'large_sigz_spectra.txt'), 'w');
fprintf(fid, '%% lambda_nm, d2W/dwdOmega [J/(rad/s)/sr] x pol. at sigma_z = %s nm; incoherent level = last column\n', mat2str(sel));
fclose(fid); writematrix([T, (pref*I/nm)'], fullfile(res_dir, 'large_sigz_spectra.txt'), 'Delimiter', ' ', 'WriteMode', 'append');
save(fullfile(res_dir, 'large_sigz_summary.mat'), 'sz', 'bands', 'bn', 'W', 'Wsd', 'Wcoh', 'Winc', 'sens', 'sens_sd', 'lam', 'Sall', 'gridchk');

% figures
col = lines(3);
fh = figure('Position', [50 50 1500 480], 'Visible', 'off');
subplot(1,3,1);
for b = [1 2]
    fill([sz sz(end:-1:1)]*1e9, [W(b,:)+Wsd(b,:), max(W(b,end:-1:1)-Wsd(b,end:-1:1), 1e-30)], col(b,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off'); hold on;
    semilogy(sz*1e9, W(b,:), '-', 'Color', col(b,:), 'DisplayName', bn{b});
    semilogy(sz([1 end])*1e9, Winc(b)*[1 1], ':', 'Color', col(b,:), 'DisplayName', [bn{b} ', incoherent']);
end
set(gca, 'YScale', 'log'); grid on; xlabel('\sigma_z [nm]'); ylabel('W in band [J/sr]'); legend('Location', 'northeast', 'FontSize', 8);
title('x-polarised energy (shaded: \pm1 jackknife sd)');
subplot(1,3,2);
for b = [1 2]
    plot(sz*1e9, 100*sens(b,:), '-', 'Color', col(b,:), 'DisplayName', bn{b}); hold on;
end
grid on; xlabel('\sigma_z [nm]'); ylabel('-d ln W / d\sigma_z [% per nm]'); legend('Location', 'northeast'); title('sensitivity');
subplot(1,3,3);
cs = parula(numel(sel) + 1);
for i = 1:numel(sel); loglog(lam*1e9, Sall(is(i),:), 'Color', cs(i,:), 'DisplayName', sprintf('%d nm', sel(i))); hold on; end
loglog(lam*1e9, pref*I/nm, 'k:', 'LineWidth', 1.5, 'DisplayName', 'incoherent');
xline(280, 'k--', 'HandleVisibility', 'off'); xline(300, 'k--', 'HandleVisibility', 'off');
grid on; xlabel('\lambda [nm]'); ylabel('d^2W/d\omegad\Omega [J/(rad/s)/sr]'); legend('Location', 'eastoutside', 'FontSize', 7);
ylim([min(pref*I/nm)/10, max(Sall(:))*3]); title('spectra (label: \sigma_z)');
sgtitle('x = 0, y = -2.0 mm, z = 5.9 m, x polariser, 100 pC (expected value over the longitudinal positions)');
savefig(fh, fullfile(fig_dir, 'large_sigz.fig')); exportgraphics(fh, fullfile(fig_dir, 'large_sigz.png'), 'Resolution', 150);
