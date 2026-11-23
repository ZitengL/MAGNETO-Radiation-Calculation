%% analyze_verify_timeshift
% Results of verify_timeshift_gpt.m: direct GPT at sigma_z = 100, 300, 500 nm vs the time shift of the 10 nm GPT run,
% same 1e4 particles. Observer x = 0, y = -2.0 mm, z = 5.9 m; x component (x polariser), undulator + dipole.
% (1) per particle: || E_direct - E_shifted || / || E_direct || over the frequencies of each band
% (2) one realisation (sum over the same particles): W = int_band pref |sum_j E_j|^2 dw, direct vs shifted
% (3) self term sum_j |E_j|^2, direct vs shifted
% (4) expected-value method (analyze_large_sigz.m) evaluated from each GPT run (10, 100, 300, 500 nm) after shifting
%     it to sigma_z = 0: W(sigma_z = 10:10:1000 nm) from runs 100/300/500 nm vs from the 10 nm run.
%     The absolute values use 1e4 particles scaled to 100 pC (fields x 5, self terms x 25, nmacro x 5).
% Output: results/verify_timeshift_gpt.txt/_summary.mat, figures/verify_timeshift_gpt.png/.fig

clear; close all;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; eps0 = 8.85418781e-12;
V = load(fullfile(res_dir, 'verify_timeshift_gpt.mat')); V = V.V;
dn = find(V.done); nch = numel(dn); fprintf('%d chunks done\n', nch);
lam = V.lambda; om = V.omega; k = om/c0; pref = eps0*c0*norm(V.obs)^2/pi;
bands = V.bands; bn = {'D1 50-280 nm', 'D2 300-1000 nm'};
mbs = @(b) lam >= bands(b,1) & lam <= bands(b,2);
Wb = @(S, b) -trapz(om(mbs(b)), S(mbs(b)));
szs = V.sz_runs(2:end); ns = numel(szs);
pid = reshape((dn(:) - 1)*2000 + (1:2000), 1, []);                 % particles of the finished chunks

fid = fopen(fullfile(res_dir, 'verify_timeshift_gpt.txt'), 'w');
fprintf(fid, '%% Direct GPT vs time shift of the 10 nm GPT run, same %d particles, x = 0, y = -2.0 mm, z = 5.9 m, x component, und + dip\n', 2000*nch);
fprintf(fid, '%% (1) per particle ||E_direct - E_shifted|| / ||E_direct||: median, 99th percentile, max\n');
for s = 1:ns
    for b = 1:2
        d = V.dpp(pid, s, b);
        fprintf(fid, '%%   sz %3.0f nm, %s: median %.2e, 99%% %.2e, max %.2e\n', szs(s)*1e9, bn{b}, median(d), prctile(d, 99), max(d));
    end
end
fprintf(fid, '%% (2) one realisation, W = int pref |sum E|^2 in band [J/sr, these %d macro-particles, not scaled]: direct, shifted, rel. diff;\n', 2000*nch);
fprintf(fid, '%%     max_w ||sum E_d|^2 - |sum E_s|^2| / |sum E_d|^2 in band\n');
Rd = zeros(ns, numel(lam)); Rs = Rd; T2 = zeros(ns, 2, 4);
pn = {'und + dip', 'undulator', 'dipole'};
for s = 1:ns
    for p = 1:3
        if p == 1; ed = squeeze(sum(V.Sd(dn,s,1,1,:) + V.Sd(dn,s,2,1,:), 1)).'; es = squeeze(sum(V.Ss(dn,s,1,1,:) + V.Ss(dn,s,2,1,:), 1)).';
        else;      ed = squeeze(sum(V.Sd(dn,s,p-1,1,:), 1)).'; es = squeeze(sum(V.Ss(dn,s,p-1,1,:), 1)).'; end
        if p == 1; Rd(s,:) = pref*abs(ed).^2; Rs(s,:) = pref*abs(es).^2; end
        for b = 1:2
            wd = Wb(pref*abs(ed).^2, b); ws = Wb(pref*abs(es).^2, b);
            m = mbs(b); mx = max(abs(abs(ed(m)).^2 - abs(es(m)).^2)./abs(ed(m)).^2);
            fprintf(fid, '%%   sz %3.0f nm, %-9s, %s: %.5e  %.5e  rel %+.2e, max_w %.2e\n', szs(s)*1e9, pn{p}, bn{b}, wd, ws, ws/wd - 1, mx);
            if p == 1; T2(s,b,:) = [wd ws ws/wd-1 mx]; end
        end
    end
end
fprintf(fid, '%% (3) self term sum |E|^2 in band: direct vs shifted, rel. diff\n');
for s = 1:ns
    id = squeeze(sum(V.Id(dn,s,1,:), 1)).'; is = squeeze(sum(V.Is(dn,s,1,:), 1)).';
    fprintf(fid, '%%   sz %3.0f nm: D1 %+.2e, D2 %+.2e\n', szs(s)*1e9, Wb(is, 1)/Wb(id, 1) - 1, Wb(is, 2)/Wb(id, 2) - 1);
end
% (4) expected-value method from each run
sz = (10:10:1000)*1e-9; f = 25/nch; nm = V.nmacro*f;
We = zeros(numel(V.sz_runs), 2, numel(sz));
for r = 1:numel(V.sz_runs)
    E0 = squeeze(sum(V.E0(dn,r,1,1,:) + V.E0(dn,r,2,1,:), 1)).'*f; I = squeeze(sum(V.I0(dn,r,1,:), 1)).'*f^2;
    for i = 1:numel(sz)
        S = pref*(exp(-(k*sz(i)).^2).*(abs(E0).^2 - I) + I/nm);
        for b = 1:2; We(r,b,i) = Wb(S, b); end
    end
end
fprintf(fid, '%% (4) expected-value method (100 pC from %d particles) from each GPT run vs from the 10 nm run:\n', 2000*nch);
fprintf(fid, '%%     max over sigma_z = 10:10:1000 nm of |W_run / W_10nm - 1|, and the value at the run''s own sigma_z\n');
for r = 2:numel(V.sz_runs)
    for b = 1:2
        q = squeeze(We(r,b,:)./We(1,b,:)) - 1; [~, i0] = min(abs(sz - V.sz_runs(r)));
        fprintf(fid, '%%   run %3.0f nm, %s: max %.2e, at %3.0f nm %+.2e (W = %.4e J/sr)\n', V.sz_runs(r)*1e9, bn{b}, max(abs(q)), sz(i0)*1e9, q(i0), We(1,b,i0));
    end
end
fprintf(fid, '%% W from the 10 nm run (1e4 particles, 100 pC) at sigma_z = 100, 300, 500 nm: D1 %s, D2 %s J/sr\n', ...
    mat2str(squeeze(We(1,1,[10 30 50]))', 4), mat2str(squeeze(We(1,2,[10 30 50]))', 4));
fclose(fid); type(fullfile(res_dir, 'verify_timeshift_gpt.txt'));
save(fullfile(res_dir, 'verify_timeshift_gpt_summary.mat'), 'szs', 'lam', 'Rd', 'Rs', 'T2', 'sz', 'We', 'bands');

% figure: realisation spectra direct vs shifted and relative difference
fh = figure('Position', [50 50 1200 450], 'Visible', 'off'); col = lines(ns);
subplot(1,2,1);
for s = 1:ns
    loglog(lam*1e9, Rd(s,:), '-', 'Color', col(s,:), 'DisplayName', sprintf('%d nm, direct GPT', round(szs(s)*1e9))); hold on;
    loglog(lam*1e9, Rs(s,:), '--', 'Color', col(s,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%d nm, time shift', round(szs(s)*1e9)));
end
xline(280, 'k:', 'HandleVisibility', 'off'); xline(300, 'k:', 'HandleVisibility', 'off');
grid on; xlabel('\lambda [nm]'); ylabel('pref |\Sigma E_x|^2 [J/(rad/s)/sr]'); legend('Location', 'southeast', 'FontSize', 8);
title(sprintf('one realisation, same %d macro-particles', 2000*nch));
subplot(1,2,2);
for s = 1:ns
    semilogy(lam*1e9, abs(Rs(s,:)./Rd(s,:) - 1), '-', 'Color', col(s,:), 'DisplayName', sprintf('%d nm', round(szs(s)*1e9))); hold on;
end
set(gca, 'XScale', 'log'); grid on; xlabel('\lambda [nm]'); ylabel('|time shift / direct - 1|'); legend('Location', 'northeast');
title('relative difference');
sgtitle('Direct GPT vs time shift of the 10 nm run: x = 0, y = -2.0 mm, z = 5.9 m, x component, undulator + dipole');
savefig(fh, fullfile(fig_dir, 'verify_timeshift_gpt.fig')); exportgraphics(fh, fullfile(fig_dir, 'verify_timeshift_gpt.png'), 'Resolution', 150);
