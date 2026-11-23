%% analyze_direct_sigz
% Results of run_direct_sigz.m: DIRECT 5e4-particle calculation (new GPT run at every sigma_z, no time shift,
% no expected value). x = 0, y = -2.0 mm, z = 5.9 m, x polariser (x component, undulator + dipole), 100 pC.
% sigma_z = 10 nm: the GPT run of run_large_sigz.m (its Er at 10 nm carries no phase factor, i.e. it is the direct sum).
% Shown as computed (user, 2026-11-06): raw S = pref |sum_j E_j|^2, nothing subtracted. Reference lines only:
%   macro-particle noise floor = pref sum_j |E_j|^2 (expected value of the noise of a 5e4-particle realisation)
%   real incoherent level      = noise floor / nmacro
% Only sigma_z with all 25 chunks done are used. Error: jackknife over the 25 chunks (fields x 25/24).
% Output: results/direct_sigz/summary.txt/.mat, results/direct_sigz/spectra.txt, figures/direct_sigz.png/.fig

clear; close all;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results'); dir_ = fullfile(res_dir, 'direct_sigz');
fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; eps0 = 8.85418781e-12;
R = load(fullfile(res_dir, 'large_sigz.mat')); R = R.R;
if R.sz_real(1) ~= R.sz_gpt; error('Er(:,1) is not the GPT sigma_z'); end
lam = R.lambda; om = R.omega; pref = eps0*c0*norm(R.obs)^2/pi; nmac = R.nmacro;
bands = [50e-9 280e-9; 300e-9 1000e-9]; bn = {'D1 50-280 nm', 'D2 300-1000 nm'};
Wb = @(S, b) -trapz(om(lam >= bands(b,1) & lam <= bands(b,2)), S(lam >= bands(b,1) & lam <= bands(b,2)));

% collect: 10 nm from run_large_sigz, the rest from run_direct_sigz
sz = R.sz_gpt; E = {squeeze(R.Er(:,1,1,1,:) + R.Er(:,1,2,1,:))}; In = {squeeze(R.I(:,3,1,:))};
f = dir(fullfile(dir_, 'direct_sz*nm.mat')); skipped = {};
for i = 1:numel(f)
    D = load(fullfile(dir_, f(i).name)); D = D.D;
    if ~all(D.done); skipped{end+1} = sprintf('%.0f nm (%d/25)', D.sigma_z*1e9, sum(D.done)); continue; end %#ok<SAGROW>
    if max(abs(D.lambda - lam)) > 0 || any(D.obs ~= R.obs); error('grid or observer differs'); end
    sz(end+1) = D.sigma_z; E{end+1} = squeeze(D.Sd(:,1,1,:) + D.Sd(:,2,1,:)); In{end+1} = squeeze(D.I(:,3,1,:)); %#ok<SAGROW>
end
[sz, o] = sort(sz); E = E(o); In = In(o); ns = numel(sz); nch = size(E{1}, 1); fj = nch/(nch-1);
S = zeros(ns, numel(lam)); Sn = S; W = zeros(2, ns); Wn = W; Wsd = W;
for s = 1:ns
    S(s,:) = pref*abs(sum(E{s}, 1)).^2; Sn(s,:) = pref*sum(In{s}, 1);
    Wj = zeros(nch, 2);
    for c = 1:nch; Sj = pref*abs(fj*sum(E{s}(setdiff(1:nch, c),:), 1)).^2; for b = 1:2; Wj(c,b) = Wb(Sj, b); end; end
    for b = 1:2
        W(b,s) = Wb(S(s,:), b); Wn(b,s) = Wb(Sn(s,:), b);
        Wsd(b,s) = sqrt((nch-1)/nch*sum((Wj(:,b) - mean(Wj(:,b))).^2));
    end
end
sens = nan(2, ns); if ns > 2; for b = 1:2; sens(b,:) = -gradient(log(W(b,:)), sz)*1e-9; end; end

fid = fopen(fullfile(dir_, 'summary.txt'), 'w');
fprintf(fid, '%% DIRECT calculation, 5e4 macro-particles (nmacro = %d), new GPT run per sigma_z; x = 0, y = -2.0 mm, z = 5.9 m, x pol., und + dip\n', nmac);
fprintf(fid, '%% W = int_band pref |sum E|^2 dw [J/sr], raw (nothing subtracted); +- = jackknife sd over %d chunks\n', nch);
fprintf(fid, '%% noise = macro-particle noise floor pref sum |E_j|^2 (expected noise of one realisation); real incoherent = noise / nmacro\n');
if ~isempty(skipped); fprintf(fid, '%% not complete, not used: %s\n', strjoin(skipped, ', ')); end
fprintf(fid, '%% sigma_z_nm | W_D1 sd noise_D1 W/noise sens_%%/nm | W_D2 sd noise_D2 W/noise sens_%%/nm\n');
for s = 1:ns
    fprintf(fid, '%5.0f', sz(s)*1e9);
    for b = 1:2; fprintf(fid, ' | %.4e %.2e %.3e %8.2f %7.2f', W(b,s), Wsd(b,s), Wn(b,s), W(b,s)/Wn(b,s), 100*sens(b,s)); end
    fprintf(fid, '\n');
end
fprintf(fid, '%% real incoherent level: D1 %.3e, D2 %.3e J/sr\n', mean(Wn(1,:))/nmac, mean(Wn(2,:))/nmac);
fclose(fid); type(fullfile(dir_, 'summary.txt'));
fid = fopen(fullfile(dir_, 'spectra.txt'), 'w');
fprintf(fid, '%% lambda_nm, raw pref|sum E|^2 [J/(rad/s)/sr] at sigma_z = %s nm, then noise floor (mean over sigma_z), real incoherent\n', mat2str(round(sz*1e9)));
fclose(fid); writematrix([lam(:)*1e9, S', mean(Sn, 1)', mean(Sn, 1)'/nmac], fullfile(dir_, 'spectra.txt'), 'Delimiter', ' ', 'WriteMode', 'append');
save(fullfile(dir_, 'summary.mat'), 'sz', 'lam', 'S', 'Sn', 'W', 'Wn', 'Wsd', 'sens', 'bands', 'nmac');

col = lines(2);
fh = figure('Position', [50 50 1500 480], 'Visible', 'off');
subplot(1,3,1);
for b = 1:2
    errorbar(sz*1e9, W(b,:), Wsd(b,:), 'o-', 'Color', col(b,:), 'MarkerSize', 4, 'DisplayName', [bn{b} ', direct']); hold on;
    plot(sz([1 end])*1e9, mean(Wn(b,:))*[1 1], '--', 'Color', col(b,:), 'DisplayName', [bn{b} ', macro-particle noise floor']);
    plot(sz([1 end])*1e9, mean(Wn(b,:))/nmac*[1 1], ':', 'Color', col(b,:), 'DisplayName', [bn{b} ', real incoherent']);
end
set(gca, 'YScale', 'log'); grid on; xlabel('\sigma_z [nm]'); ylabel('W in band [J/sr]'); legend('Location', 'northeast', 'FontSize', 7);
ylim([1e-8 1e5]);   % room for the legend above the data
title('raw energy (error bars: jackknife sd)');
subplot(1,3,2);
for b = 1:2; plot(sz*1e9, 100*sens(b,:), 'o-', 'Color', col(b,:), 'MarkerSize', 4, 'DisplayName', bn{b}); hold on; end
grid on; xlabel('\sigma_z [nm]'); ylabel('-d ln W / d\sigma_z [% per nm]'); legend('Location', 'northeast'); title('sensitivity (raw W)');
subplot(1,3,3); cs = parula(ns + 1);
for s = 1:ns; loglog(lam*1e9, S(s,:), 'Color', cs(s,:), 'DisplayName', sprintf('%.0f nm', sz(s)*1e9)); hold on; end
loglog(lam*1e9, mean(Sn, 1), 'k--', 'LineWidth', 1.5, 'DisplayName', 'noise floor');
loglog(lam*1e9, mean(Sn, 1)/nmac, 'k:', 'LineWidth', 1.5, 'DisplayName', 'real incoherent');
xline(280, 'k-.', 'HandleVisibility', 'off'); xline(300, 'k-.', 'HandleVisibility', 'off');
grid on; xlabel('\lambda [nm]'); ylabel('d^2W/d\omegad\Omega [J/(rad/s)/sr]'); legend('Location', 'eastoutside', 'FontSize', 7);
title('raw spectra (label: \sigma_z)');
sgtitle('Direct 5e4-particle calculation: x = 0, y = -2.0 mm, z = 5.9 m, x polariser, 100 pC');
savefig(fh, fullfile(fig_dir, 'direct_sigz.fig')); exportgraphics(fh, fullfile(fig_dir, 'direct_sigz.png'), 'Resolution', 150);
