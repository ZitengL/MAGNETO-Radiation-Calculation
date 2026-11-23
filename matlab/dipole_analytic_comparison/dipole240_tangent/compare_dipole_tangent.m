%% compare_dipole_tangent
% Single-electron dipole spectrum, 240 MeV, observed on the tangent to the arc at its midpoint
% (direction = half the bend angle, 40.555 mrad upwards), three calculations:
%   (1) MAGNETO, hard-edge dipole: gpt_sims/single_electron/E240MeV_dipole_L10cm_hardedge (b1 = 1e5 1/m)
%   (2) MAGNETO, Enge-fringe dipole (old model): gpt_sims/single_electron/E240MeV_dipole_L10cm (gap 1 cm)
%   (3) analytic, Jackson (14.79) / Wiedemann (25.88) at theta = 0 (in the bend plane):
%       d2W/dw dOmega = 3 e^2 gamma^2/(16 pi^3 eps0 c) (w/wc)^2 K_{2/3}^2(w/(2 wc)),
%       wc = 3 gamma^3 c/(2 rho), rho = gamma beta m c/(e B0), B0 = 0.65 T (uniform circular arc, far field)
% Both GPT dipoles: centre z = 1 m, effective length 0.1 m, 0.65 T, same field integral and bend angle;
% the only difference is the fringe field.
% Observers on the midpoint tangent of the hard-edge trajectory (the Enge tangent is 33 um away at 11 m):
%   '11m' : intersection with the plane z = 11 m (10.008 m from the tangent point)
%   '1km' : 1000 m from the tangent point along the tangent (far field check)
% Normalisation: per steradian seen from the tangent point P; MAGNETO x (|obs - P|/|obs|)^2.
% Checks: tout (3.75 vs 1.875 um, both dipoles), hard edge (b1 = 1e5 vs 1e6), polarisation (pi = x ~ 0).
% 15-150 nm, 5000 log points.
% Outputs: results/compare_dipole_tangent.mat, results/dipole_tangent_<obs>.txt,
%          figures/compare_dipole_tangent_<obs>.(png|fig)

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
if ~exist(res_dir, 'dir'); mkdir(res_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
c0 = 299792458; eps0 = 8.85418781e-12; me_g = 9.10953e-31; qe_g = 1.6021892e-19;   % GPT constants for e, m
B0 = 0.65;
lam_nm = logspace(log10(15), log10(150), 5000); lam = lam_nm*1e-9; w = 2*pi*c0 ./ lam;

sdir = fullfile(root_dir, 'gpt_sims', 'single_electron');
hd = fullfile(sdir, 'E240MeV_dipole_L10cm_hardedge'); en = fullfile(sdir, 'E240MeV_dipole_L10cm');
F = struct('hard',      fullfile(hd, 'result.gdf'), ...
           'hard_b1e6', fullfile(hd, 'b1_check', 'result_b1e6.gdf'), ...
           'hard_dz1p875', fullfile(hd, 'tout_convergence', 'result_dz1p875um.gdf'), ...
           'enge',      fullfile(en, 'result.gdf'), ...
           'enge_dz1p875', fullfile(en, 'tout_convergence', 'result_dz1p875um.gdf'));
fn = fieldnames(F);
for i = 1:numel(fn); T.(fn{i}) = field_region(read_gpt_trajectory(F.(fn{i}))); end

%% midpoint tangent of the hard-edge trajectory
tr = T.hard; ang = atan2(tr.b(2,:), tr.b(3,:)); th = ang(end); k = find(ang >= th/2, 1);
s = (th/2 - ang(k-1)) / (ang(k) - ang(k-1)); P = (tr.r(:, k-1) + s*(tr.r(:, k) - tr.r(:, k-1)))';
d = [0 sin(th/2) cos(th/2)];
obs.o11m = [0, P(2) + (11 - P(3))*tan(th/2), 11];
obs.o1km = P + 1000*d;
gam = tr.G(1); GB = sqrt(gam^2 - 1);
rho = GB*me_g*c0/(qe_g*B0); wc = 3*gam^3*c0/(2*rho);
fprintf('bend %.6f mrad, tangent point P = (0, %.6f mm, %.6f m), direction %.6f mrad\n', th*1e3, P(2)*1e3, P(3), th/2*1e3);
fprintf('gamma %.8f, rho %.6f m, lambda_c %.4f nm; observers: 11 m -> (0, %.6f, 11) m; 1 km -> (0, %.4f, %.4f) m\n', ...
    gam, rho, 2*pi*c0/wc*1e9, obs.o11m(2), obs.o1km(2), obs.o1km(3));

%% analytic
yv = w / wc;
S_an = 3*qe_g^2*gam^2/(16*pi^3*eps0*c0) * yv.^2 .* besselk(2/3, yv/2).^2;

%% MAGNETO
on = fieldnames(obs); out = struct();
for q = 1:numel(on)
    O = obs.(on{q}); fac = (norm(O - P)/norm(O))^2;
    for i = 1:numel(fn)
        sp = emission_time_spectrum(T.(fn{i}), O, lam, 'fB');
        S.(fn{i}) = fac * sp.d2W_dwdOmega.total;
        px.(fn{i}) = max(sp.d2W_dwdOmega.x ./ sp.d2W_dwdOmega.total);
    end
    r_h = S.hard ./ S_an - 1; r_e = S.enge ./ S_an - 1; r_eh = S.enge ./ S.hard - 1;
    fprintf('\n=== %s: R_src = %.4f m\n', on{q}, norm(O - P));
    fprintf('checks: |hard 3.75/1.875 um - 1| max %.1e; |Enge 3.75/1.875 um - 1| max %.1e; |hard b1 1e5/1e6 - 1| max %.1e; pi (x) / total max: hard %.1e, Enge %.1e\n', ...
        max(abs(S.hard./S.hard_dz1p875 - 1)), max(abs(S.enge./S.enge_dz1p875 - 1)), max(abs(S.hard./S.hard_b1e6 - 1)), px.hard, px.enge);
    for band = [15 25; 25 40; 40 70; 70 150]'
        sel = lam_nm >= band(1) & lam_nm <= band(2);
        fprintf('  %3d-%3d nm: hard/analytic - 1 in [%+.2e, %+.2e]; Enge/analytic - 1 in [%+.2e, %+.2e]; Enge/hard - 1 in [%+.2e, %+.2e]\n', band, ...
            min(r_h(sel)), max(r_h(sel)), min(r_e(sel)), max(r_e(sel)), min(r_eh(sel)), max(r_eh(sel)));
    end
    out.(on{q}) = struct('obs', O, 'R_src', norm(O - P), 'S_analytic', S_an, 'S_MAGNETO_hardedge', S.hard, 'S_MAGNETO_enge', S.enge, ...
        'S_hard_dz1p875', S.hard_dz1p875, 'S_hard_b1e6', S.hard_b1e6, 'S_enge_dz1p875', S.enge_dz1p875, ...
        'ratio_hard_analytic_m1', r_h, 'ratio_enge_analytic_m1', r_e, 'ratio_enge_hard_m1', r_eh);

    fid = fopen(fullfile(res_dir, ['dipole_tangent_' on{q}(2:end) '.txt']), 'wt');
    fprintf(fid, '%% 240 MeV single electron, dipole 0.1 m / 0.65 T, observer on the arc-midpoint tangent at (%.6g, %.6g, %.6g) m (%s)\n', O, on{q}(2:end));
    fprintf(fid, '%% d2W/domega dOmega [J/(rad/s)/sr] per sr seen from the tangent point (R_src = %.4f m), total polarisation\n', norm(O - P));
    fprintf(fid, '%% Columns: 1 lambda [nm]  2 omega [rad/s]  3 analytic Jackson (14.79) theta = 0  4 MAGNETO hard-edge dipole  5 MAGNETO Enge-fringe dipole (gap 1 cm)\n');
    fprintf(fid, '%-14.8f %-14.6e %-14.6e %-14.6e %-14.6e\n', [lam_nm; w; S_an; S.hard; S.enge]);
    fclose(fid);

    f = figure('Position', [100 100 1000 750]);
    subplot(3,1,[1 2]);
    loglog(lam_nm, S_an, 'k-', 'LineWidth', 1.6); hold on;
    loglog(lam_nm, S.hard, 'b--', 'LineWidth', 1.2); loglog(lam_nm, S.enge, 'r:', 'LineWidth', 1.4);
    grid on; xlim([15 150]); set(gca, 'XTick', [15 20 30 40 50 70 100 150]);
    legend('analytic, Jackson (14.79), \theta = 0', 'MAGNETO, hard-edge dipole', 'MAGNETO, Enge fringe (gap 1 cm)', 'Location', 'southeast');
    ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]');
    title(sprintf('240 MeV single electron, dipole 0.1 m / 0.65 T (\\lambda_c = %.1f nm): arc-midpoint tangent, %s', 2*pi*c0/wc*1e9, ...
        strrep(strrep(on{q}(2:end), 'm', ' m'), 'k m', ' km')));
    subplot(3,1,3);
    semilogx(lam_nm, r_h, 'b-', lam_nm, r_e, 'r-', lam_nm, r_eh, 'm-'); grid on; xlim([15 150]);
    set(gca, 'XTick', [15 20 30 40 50 70 100 150]);
    legend('hard edge / analytic - 1', 'Enge / analytic - 1', 'Enge / hard edge - 1', 'Location', 'best');
    xlabel('\lambda [nm]'); ylabel('relative difference');
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(f, fullfile(fig_dir, ['compare_dipole_tangent_' on{q}(2:end) '.png']), 'Resolution', 200);
    savefig(f, fullfile(fig_dir, ['compare_dipole_tangent_' on{q}(2:end) '.fig'])); close(f);
end
out.lambda_nm = lam_nm; out.omega = w; out.tangent_point = P; out.tangent_dir = d; out.bend_rad = th;
out.gamma = gam; out.rho = rho; out.omega_c = wc; out.B0 = B0;
out.note = 'per sr seen from the tangent point; analytic Jackson (14.79) theta = 0; hard edge b1 = 1e5 1/m, Enge b1 = 200 1/m';
save(fullfile(res_dir, 'compare_dipole_tangent.mat'), 'out');

function traj = field_region(traj)
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
