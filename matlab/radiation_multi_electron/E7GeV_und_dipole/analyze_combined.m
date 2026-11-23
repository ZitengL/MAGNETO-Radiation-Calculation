%% analyze_combined
% Results of run_combined.m (undulator turned over + dipole bent up, 5e4 macro-particles, 100 pC), observers
% x = 0, y = -1 ... -3 mm, z = 5.9 m, 50-300 nm, sigma_z = 10:5:50 nm (time shift from the 10 nm GPT run).
%  - spectra and band energies per part (undulator / dipole / sum) and polarisation (x = what an x polariser passes, y)
%  - dipole contamination of the x-polarised signal: (W_tot,x - W_und,x) / W_und,x (includes the interference)
%  - polariser requirement: extinction ratio eps so that eps * W_tot,y = 1 % of W_tot,x
%  - self terms (macro-particle noise) and jackknife errors over the chunks (fields x nch/(nch-1))
%  - analytic positions-only model (3D Gaussian form factor, theta from the undulator centre, sigma = 20 um),
%    real N: N S1 + N(N-1) |F|^2 S1, S1 = single ideal electron (und + dip, x pol.) at the same observer
%  - check 1: undulator part at y = -1 mm vs the direct GPT sigma_z sweep (E7GeV_gap6und_sigz_sweep, rotated image)
%  - check 2: chunks 1-3 at sigma_z = 50 nm, time shift vs direct GPT run (direct50_combined.mat)
% Output: results/combined_summary.txt/.mat, results/combined_spectra_y*.txt, figures/combined_*.png

clear; close all;
this_dir = fileparts(mfilename('fullpath')); root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; eps0 = 8.85418781e-12; N = 100e-12/1.602176634e-19;
R = load(fullfile(res_dir, 'combined.mat')); R = R.R;
dn = find(R.done); nch = numel(dn); nall = numel(R.done);
if nch < nall; fprintf('WARNING: only %d of %d chunks done; fields scaled by %d/%d\n', nch, nall, nall, nch); end
sc = nall/nch;                                                     % keep the total charge when chunks are missing
om = R.omega; lam = R.lambda; szl = R.sigma_z; ns = numel(szl); no = size(R.obs, 1); ys = R.obs(:,2)';
W = @(S) -trapz(om, S, 2);
Ec = R.Ec(dn,:,:,:,:,:)*sc; Ic = R.Ic(dn,:,:,:,:)*sc^2/sc;          % self terms scale with the number of particles (x sc)
Eall = squeeze(sum(Ec, 1));                                        % sz x obs x part x comp x omega
Iall = squeeze(sum(Ic, 1));                                        % obs x part x comp x omega

% single ideal electron (S1 for the analytic model)
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole', 'ideal_electron');
tr = read_gpt_trajectory(fullfile(gdir, 'result_sm_dz60.gdf'));
S1 = zeros(no, numel(lam));
for o = 1:no
    s = emission_time_spectrum(tr, R.obs(o,:), lam, 'fB', norm(R.obs(o,:))/c0);
    S1(o,:) = eps0*c0*norm(R.obs(o,:))^2/pi*abs(s.Ew(1,:)).^2;      % x polarisation, und + dip
end

T = []; k = om/c0;
for o = 1:no
    pref = eps0*c0*norm(R.obs(o,:))^2/pi;
    th = R.obs(o,1:2)/(R.obs(o,3) - 0.4);
    Fpos = exp(-k.^2*(20e-6)^2*sum(th.^2));
    Iself_x = pref*squeeze(Iall(o,3,1,:))';                        % self terms of the x-polarised sum (sz independent)
    for s = 1:ns
        Eu = squeeze(Eall(s,o,1,:,:)); Ed = squeeze(Eall(s,o,2,:,:)); Et = Eu + Ed;
        Sux = pref*abs(Eu(1,:)).^2; Sdx = pref*abs(Ed(1,:)).^2; Stx = pref*abs(Et(1,:)).^2; Sty = pref*abs(Et(2,:)).^2;
        Suy = pref*abs(Eu(2,:)).^2; Sdy = pref*abs(Ed(2,:)).^2;
        % jackknife of W_tot,x
        jk = zeros(1, nch);
        for j = 1:nch
            Ej = (squeeze(sum(Ec(:,s,o,1,1,:) + Ec(:,s,o,2,1,:), 1)) - squeeze(Ec(j,s,o,1,1,:) + Ec(j,s,o,2,1,:)))*nch/(nch-1);
            jk(j) = W(pref*abs(Ej(:)').^2);
        end
        jsd = sqrt((nch-1)/nch*sum((jk - mean(jk)).^2));
        Fz = exp(-k.^2*szl(s)^2);
        WA = W(N*S1(o,:) + N*(N-1)*Fpos.*Fz.*S1(o,:));
        T(end+1,:) = [1e3*ys(o), 1e9*szl(s), W(Stx), jsd, W(Iself_x), W(Sux), W(Sdx), (W(Stx) - W(Sux))/W(Sux), W(Sty), W(Suy), W(Sdy), ...
                      0.01*W(Stx)/W(Sty), WA]; %#ok<SAGROW>
        Spec{o, s} = [Sux; Sdx; Stx; Sty]; %#ok<SAGROW>
    end
end
hdr = {'y_mm', 'sigma_z_nm', 'W_tot_x', 'W_tot_x_jackknife_sd', 'W_self_x', 'W_und_x', 'W_dip_x', 'dipole_contamination_x', ...
       'W_tot_y', 'W_und_y', 'W_dip_y', 'polariser_extinction_for_1pct', 'W_analytic_positions_only_x'};
fid = fopen(fullfile(res_dir, 'combined_summary.txt'), 'wt');
fprintf(fid, '%% undulator (turned over) + dipole (bent up), 100 pC, 5e4 macro-particles (%d of %d chunks), x = 0, z = 5.9 m, energy in 50-300 nm [J/sr]\n', nch, nall);
fprintf(fid, '%% x = what an x polariser passes. dipole_contamination_x = (W_tot_x - W_und_x)/W_und_x (incl. interference);\n');
fprintf(fid, '%% polariser_extinction_for_1pct = extinction ratio at which the leaked y light is 1 %% of W_tot_x; W_self_x = macro-particle self terms (noise level)\n');
fprintf(fid, '%% %s\n', strjoin(hdr, ' ')); fprintf(fid, [repmat('%.6e ', 1, numel(hdr)) '\n'], T'); fclose(fid);
for o = 1:no
    M = lam(:)*1e9;
    for s = 1:ns; M = [M, Spec{o,s}']; end %#ok<AGROW>
    fid = fopen(fullfile(res_dir, sprintf('combined_spectra_y%.1fmm.txt', 1e3*ys(o))), 'wt');
    fprintf(fid, '%% y = %.1f mm, d2W/dwdOmega [J/(rad/s)/sr]; for each sigma_z = 10:5:50 nm four columns: und_x dip_x tot_x tot_y\n', 1e3*ys(o));
    fprintf(fid, [repmat('%.6e ', 1, size(M,2)) '\n'], M'); fclose(fid);
end
save(fullfile(res_dir, 'combined_summary.mat'), 'T', 'hdr', 'Spec', 'lam', 'szl', 'ys', 'S1', 'nch');

% print
fprintf('%s\n', strjoin(hdr, ' | '));
for i = 1:size(T, 1); fprintf('%5.1f %3.0f | %.4e +- %.1e | self %.2e | und %.4e dip %.2e | contam %+.2e | y: %.3e (und %.3e dip %.3e) | eps %.1e | analytic %.4e\n', T(i,:)); end

% check 1: undulator part at y = -1 mm vs direct sigma_z sweep (rotated image: (Ex,Ey,Ez) -> (-Ex,-Ey,Ez))
o1 = find(abs(ys + 1e-3) < 1e-9);
for s = 1:ns
    f = fullfile(this_dir, '..', 'E7GeV_gap6und_sigz_sweep', 'results', sprintf('sigz%02.0fnm.mat', 1e9*szl(s)));
    if ~exist(f, 'file'); continue; end
    Q = load(f); Q = Q.R; dd = 0; mm = 0;
    for j = dn
        En = squeeze(R.Ec(j,s,o1,1,:,:)); Eo = diag([-1 -1 1])*squeeze(Q.Ec(j,:,:));
        dd = max(dd, max(abs(En - Eo), [], 'all')); mm = max(mm, max(abs(Eo), [], 'all'));
    end
    fprintf('check 1, sigma_z %2.0f nm: undulator part (time shift, dz 120 um) vs direct sweep (dz 60 um), per chunk: max |dE| / max |E| = %.2e\n', 1e9*szl(s), dd/mm);
end
% check 2: direct GPT at 50 nm
f = fullfile(res_dir, 'direct50_combined.mat');
if exist(f, 'file')
    D = load(f); D = D.R; s50 = find(abs(szl - 50e-9) < 1e-12);
    for j = find(D.done)
        for part = 1:2
            En = squeeze(R.Ec(j,s50,:,part,:,:)); Ed = squeeze(D.Ec(j,1,:,part,:,:));
            fprintf('check 2, chunk %d, part %d: time shift vs direct GPT at 50 nm: max |dE| / max |E| = %.2e\n', j, part, max(abs(En - Ed), [], 'all')/max(abs(Ed), [], 'all'));
        end
    end
end

% figures
cm = lines(no);
f = figure('Position', [50 50 1100 420]);
subplot(1,2,1); hold on;
for o = 1:no
    r = T(:,1) == 1e3*ys(o);
    errorbar(T(r,2), T(r,3), T(r,4), 'o-', 'Color', cm(o,:), 'MarkerFaceColor', cm(o,:), 'DisplayName', sprintf('y = %.1f mm, MAGNETO', 1e3*ys(o)));
    plot(T(r,2), T(r,13), '--', 'Color', cm(o,:), 'HandleVisibility', 'off');
end
set(gca, 'YScale', 'log'); grid on; xlabel('\sigma_z [nm]'); ylabel('x-polarised energy, 50-300 nm [J/sr]'); legend('Location', 'best');
title('undulator + dipole, x polariser (dashed: analytic, positions only)');
subplot(1,2,2); hold on;
for o = 1:no
    r = T(:,1) == 1e3*ys(o);
    plot(T(r,2), T(r,3)/T(find(r,1),3), 'o-', 'Color', cm(o,:), 'MarkerFaceColor', cm(o,:), 'DisplayName', sprintf('y = %.1f mm', 1e3*ys(o)));
end
set(gca, 'YScale', 'log'); grid on; xlabel('\sigma_z [nm]'); ylabel('W_x / W_x(\sigma_z = 10 nm)'); legend('Location', 'best');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'combined_W_vs_sigz.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'combined_W_vs_sigz.fig')); close(f);

f = figure('Position', [50 50 1400 700]);
for o = 1:no
    subplot(2,3,o); hold on; X = Spec{o,1};
    plot(lam*1e9, X(1,:), 'b-', 'DisplayName', 'undulator, x'); plot(lam*1e9, X(2,:), 'r-', 'DisplayName', 'dipole, x');
    plot(lam*1e9, X(3,:), 'k-', 'DisplayName', 'sum, x'); plot(lam*1e9, X(4,:), 'm-', 'DisplayName', 'sum, y (blocked)');
    pref = eps0*c0*norm(R.obs(o,:))^2/pi; plot(lam*1e9, pref*squeeze(Iall(o,3,1,:))', 'k:', 'DisplayName', 'self terms, x');
    set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; xlim([50 300]); xlabel('\lambda [nm]'); ylabel('d^2W/d\omegad\Omega [J/(rad/s)/sr]');
    title(sprintf('\\sigma_z = 10 nm, x = 0, y = %.1f mm', 1e3*ys(o))); if o == 1; legend('Location', 'best'); end
end
subplot(2,3,6); hold on; r = T(:,2) == 10;
plot(T(r,1), abs(T(r,8)), 'ko-', 'DisplayName', '|dipole contamination|, x');
plot(T(r,1), T(r,12), 'mo-', 'DisplayName', 'polariser extinction for 1 % leak');
plot(T(r,1), T(r,5)./T(r,3), 'bo-', 'DisplayName', 'self terms / W_x');
set(gca, 'YScale', 'log'); grid on; xlabel('y [mm]'); title('\sigma_z = 10 nm'); legend('Location', 'best');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'combined_spectra.png'), 'Resolution', 150); savefig(f, fullfile(fig_dir, 'combined_spectra.fig')); close(f);
