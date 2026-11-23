%% compare_undulator
% Single-electron undulator spectra, 240 MeV: analytic (Wiedemann) vs MAGNETO.
%   analytic : (26.46) full double Bessel sum (main), (26.51) per harmonic (check near the peaks),
%              (26.59)-(26.61) on-axis peak values; N = 15, K = 0.980309, lambda_u = 30.472767 mm,
%              gamma of the GPT run; third integral (26.36) not included (user's decision)
%   MAGNETO  : emission_time_spectrum.m (fB) on
%              (a) ideal undulator: gpt_sims/single_electron/E240MeV_unduplan_ideal (unduplan, Nu = 16 with
%                  1/4, 3/4 end half periods = 15 effective periods)
%              (b) real Radia map: gpt_sims/single_electron/E240MeV_und15p_paraxial
% Observers (wiggle plane, phi = 0), angle theta_u seen from the undulator centre (z = 1 m):
%   '11m_theta0'    (0, 0, 11) m                                    theta_u = 0
%   '11m_theta0p5'  11 m (sin 0.5 mrad, 0, cos 0.5 mrad)              theta_u = 0.5500 mrad
%   '1km_theta0'    (0, 0, 1) + 1000 m (0, 0, 1)                     theta_u = 0         (far field)
%   '1km_theta0p5'  (0, 0, 1) + 1000 m (sin th_u, 0, cos th_u)       theta_u = 0.5500 mrad (far field)
% Normalisation: the analytic formula is per steradian seen from the source; MAGNETO is per steradian
% of the sphere |obs| centred at the origin, so MAGNETO is multiplied by (R_src/|obs|)^2 with
% R_src = |obs - (0,0,1 m)|. 15-150 nm, 5000 log points.
% Outputs: results/compare_undulator.mat, results/undulator_<case>.txt,
%          figures/compare_<case>.(png|fig) (spectra + difference)

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
if ~exist(res_dir, 'dir'); mkdir(res_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
c0 = 299792458;
K = 0.980309; lamu = 0.030472767; Np = 15; zu = 1.0;
lam_nm = logspace(log10(15), log10(150), 5000); lam = lam_nm*1e-9;
nw = round(0.01 / mean(diff(log(lam_nm)))); env = @(S) movmax(S, 2*nw + 1);

sdir = fullfile(root_dir, 'gpt_sims', 'single_electron');
T.ideal = field_region(read_gpt_trajectory(fullfile(sdir, 'E240MeV_unduplan_ideal', 'result.gdf')));
T.real  = field_region(read_gpt_trajectory(fullfile(sdir, 'E240MeV_und15p_paraxial', 'result.gdf')));
gam = T.ideal.G(1);

thu5 = atan(11*sin(0.5e-3) / (11*cos(0.5e-3) - 1));
cases = struct('name', {'11m_theta0', '11m_theta0p5', '1km_theta0', '1km_theta0p5'}, ...
    'obs', {[0 0 11], 11*[sin(0.5e-3) 0 cos(0.5e-3)], [0 0 zu+1000], [0 0 zu] + 1000*[sin(thu5) 0 cos(thu5)]}, ...
    'thu', {0, thu5, 0, thu5}, 'label', {'11 m, \theta_x = 0', '11 m, \theta_x = 0.5 mrad (\theta_u = 0.550 mrad)', ...
    '1 km (far field), \theta_u = 0', '1 km (far field), \theta_u = 0.550 mrad'});

out = struct();
for q = 1:numel(cases)
    cs = cases(q);
    a = undulator_analytic(lam, gam, K, Np, lamu, cs.thu, 0);
    Rsrc = norm(cs.obs - [0 0 zu]); conv = (Rsrc / norm(cs.obs))^2;
    si = emission_time_spectrum(T.ideal, cs.obs, lam, 'fB');
    sr = emission_time_spectrum(T.real,  cs.obs, lam, 'fB');
    Si = conv * si.d2W_dwdOmega.total; Sr = conv * sr.d2W_dwdOmega.total;
    Sa = a.W46.total; S51 = a.W51.total;
    e_i = abs(Si - Sa) ./ env(Sa); e_r = abs(Sr - Sa) ./ env(Sa); e_51 = abs(S51 - Sa) ./ env(Sa);
    fprintf('\n=== %s: R_src = %.4f m, theta_u = %.4f mrad; analytic lambda_1 = %.4f nm; Bessel |n| <= %d\n', cs.name, Rsrc, cs.thu*1e3, a.lambda1*1e9, a.nmax);
    fprintf('MAGNETO polarisation: ideal max y/total %.1e, z/total %.1e; real max y/total %.1e\n', ...
        max(si.d2W_dwdOmega.y ./ si.d2W_dwdOmega.total), max(si.d2W_dwdOmega.z ./ si.d2W_dwdOmega.total), max(sr.d2W_dwdOmega.y ./ sr.d2W_dwdOmega.total));
    fprintf('|diff|/envelope of the analytic (max / median): ideal %.2e / %.1e, real %.2e / %.1e, (26.51) vs (26.46) %.2e / %.1e\n', ...
        max(e_i), median(e_i), max(e_r), median(e_r), max(e_51), median(e_51));
    if cs.thu == 0; kl = [1 3 5]; else; kl = 1:5; end
    lam_k = a.lambda1*1e9 ./ kl;
    fprintf('  k  lambda_k [nm] | peak position [nm]: (26.46)  ideal  real | peak value / (26.46) peak: ideal  real  (26.51)  (26.61) | ideal peak/max\n');
    pk = zeros(numel(kl), 9);
    for j = 1:numel(kl)
        win = abs(lam_nm - lam_k(j)) < 0.04*lam_k(j);
        [la, va] = peak_refined(lam_nm, Sa, win); [li, vi] = peak_refined(lam_nm, Si, win);
        [lr, vr] = peak_refined(lam_nm, Sr, win); [~, v51] = peak_refined(lam_nm, S51, win);
        p61 = NaN; if cs.thu == 0; p61 = a.peak61(a.k_odd == kl(j)); end
        pk(j, :) = [kl(j) lam_k(j) la li lr vi/va vr/va v51/va p61/va];
        fprintf('  %d  %9.4f   | %9.4f %9.4f %9.4f | %7.4f %7.4f %7.4f %7.4f | %.1e\n', pk(j, 1:9), vi/max(Si));
    end
    sel = Sa > 1e-2*max(Sa);                                % where the spectrum is not tiny
    fprintf('  where S_analytic > 1%% of its max: ideal/analytic - 1 in [%+.3f, %+.3f], real/analytic - 1 in [%+.3f, %+.3f]\n', ...
        min(Si(sel)./Sa(sel) - 1), max(Si(sel)./Sa(sel) - 1), min(Sr(sel)./Sa(sel) - 1), max(Sr(sel)./Sa(sel) - 1));

    r = struct('obs', cs.obs, 'theta_u', cs.thu, 'R_src', Rsrc, 'S_analytic_2646', Sa, 'S_analytic_2651', S51, ...
        'S_analytic_x', a.W46.x, 'S_MAGNETO_ideal', Si, 'S_MAGNETO_real', Sr, 'err_env_ideal', e_i, 'err_env_real', e_r, ...
        'peaks', pk, 'peak61', a.peak61, 'lambda1_nm', a.lambda1*1e9);
    r.peaks_columns = 'k, lambda_k [nm], peak position (26.46), ideal, real [nm], peak value / (26.46) peak: ideal, real, (26.51), (26.61)';
    out.(['c_' cs.name]) = r;

    fid = fopen(fullfile(res_dir, ['undulator_' cs.name '.txt']), 'wt');
    fprintf(fid, '%% Single electron, 240 MeV, 15-period undulator (K = %.6f, lambda_u = %.6f mm), observer %s: obs = (%.6g, %.6g, %.6g) m, theta_u = %.4f mrad\n', K, lamu*1e3, cs.name, cs.obs, cs.thu*1e3);
    fprintf(fid, '%% d2W/domega dOmega [J/(rad/s)/sr] per steradian seen from the undulator centre (R_src = %.4f m), total polarisation\n', Rsrc);
    fprintf(fid, '%% Columns: 1 lambda [nm]  2 omega [rad/s]  3 analytic (26.46)  4 analytic (26.51), valid near peaks only  5 MAGNETO ideal (unduplan)  6 MAGNETO real map\n');
    fprintf(fid, '%-14.8f %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e\n', [lam_nm; 2*pi*c0./lam; Sa; S51; Si; Sr]);
    fclose(fid);

    f = figure('Position', [100 100 1000 750]);
    subplot(3,1,[1 2]);
    loglog(lam_nm, Sa, 'k-', 'LineWidth', 1.6); hold on;
    loglog(lam_nm, Si, 'b--', 'LineWidth', 1.1); loglog(lam_nm, Sr, 'r:', 'LineWidth', 1.2);
    for j = 1:numel(kl); xline(lam_k(j), 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); end
    grid on; xlim([15 150]); set(gca, 'XTick', [15 20 30 40 50 70 100 150]); ylim(max(Sa)*[1e-6 3]);
    legend('analytic, Wiedemann (26.46), N = 15', 'MAGNETO, ideal field (unduplan)', 'MAGNETO, real map (Radia)', 'Location', 'southeast');
    ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]'); title(['240 MeV single electron, 15-period undulator: ' cs.label]);
    subplot(3,1,3);
    semilogy(lam_nm, e_i, 'b-', lam_nm, e_r, 'r-'); set(gca, 'XScale', 'log'); grid on; xlim([15 150]);
    set(gca, 'XTick', [15 20 30 40 50 70 100 150]); ylim([1e-5 10]);
    legend('|ideal - analytic|', '|real - analytic|', 'Location', 'best');
    xlabel('\lambda [nm]'); ylabel('diff / local envelope');
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(f, fullfile(fig_dir, ['compare_' cs.name '.png']), 'Resolution', 200); savefig(f, fullfile(fig_dir, ['compare_' cs.name '.fig']));
    close(f);
end
out.lambda_nm = lam_nm; out.omega = 2*pi*c0./lam; out.params = struct('K', K, 'lambda_u', lamu, 'N_analytic', Np, 'gamma', gam, ...
    'note', 'S per sr seen from the undulator centre; analytic without the third integral (26.36)');
save(fullfile(res_dir, 'compare_undulator.mat'), 'out');

function [lp, Sp] = peak_refined(lam, S, win)
[~, j] = max(S .* win);
x = log(lam(j-1:j+1)); y = log(S(j-1:j+1)); p = polyfit(x - x(2), y, 2);
dx = -p(2)/(2*p(1)); lp = exp(x(2) + dx); Sp = exp(polyval(p, dx));
end

function traj = field_region(traj)
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
