%% compare_dipole_onaxis
% 240 MeV short dipole (rectmagnet L = 0.1 m, B0 = 0.65 T, Enge gap 1 cm, centre z = 1 m):
% MAGNETO on-axis spectrum at (0, 0, 11) m (matlab/radiation_dipole240_single_electron)
% versus the analytic angle-resolved bending-magnet spectrum, Jackson (14.79) /
% Wiedemann (25.88), in the orbit plane (theta = 0, sigma polarisation only):
%   d2W/dw dOmega = 3 e^2 gamma^2 / (16 pi^3 eps0 c) * y^2 * K_{2/3}(y/2)^2,  y = w / w_c,
%   w_c = 3 gamma^3 c / (2 rho),  rho = gamma beta m c / (e B0)
% The formula assumes a long circular arc in a uniform field around the tangent point.
% The on-axis observer looks along the incoming straight line: the electron's direction
% points at the observer only at the very start of the bend (in the entrance fringe
% field), with no arc before it, so a large deviation is expected (edge effect). This
% script quantifies it.
% Normalisation: MAGNETO gives per sr of the R = 11 m sphere about the origin; the
% analytic formula is per sr seen from the source. MAGNETO values are converted with
% (R_src / 11 m)^2, R_src = distance from the tangent point to the observer.
% Outputs: results/compare_dipole_onaxis.(mat|txt), figures/compare_dipole_onaxis.(png|fig)

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
if ~exist(res_dir, 'dir'); mkdir(res_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
e0 = 1.60217663e-19; c0 = 299792458; eps0 = 8.85418781e-12; me = 9.1093837e-31;

S = load(fullfile(root_dir, 'matlab', 'radiation_dipole240_single_electron', 'results', 'dipole240_onaxis.mat')); d = S.d;
lam = d.lambda_nm; w = d.omega;
Sm_tot = d.d2W_dwdOmega; Sm_xyz = d.d2W_dwdOmega_xyz;
obs = [0 0 11];

traj = read_gpt_trajectory(fullfile(root_dir, 'gpt_sims', 'single_electron', 'E240MeV_dipole_L10cm', 'result.gdf'));
z = traj.r(3,:); y = traj.r(2,:); B = traj.fB(1,:);
gam = traj.G(1); GB = sqrt(gam^2 - 1);
[~, ic] = min(abs(z - 1.0)); B0 = abs(B(ic));
rho = GB*me*c0/(e0*B0); wc = 3*gam^3*c0/(2*rho);
fprintf('gamma %.4f, B0 %.5f T, rho %.4f m, lambda_c %.3f nm, bend angle %.3f mrad\n', gam, B0, rho, 2*pi*c0/wc*1e9, ...
    1e3*atan2(traj.b(2,end), traj.b(3,end)));
% tangent point: where the electron direction points at the observer
alpha = atan2(traj.b(2,:), traj.b(3,:));            % direction of motion
phi   = atan2(obs(2) - y, obs(3) - z);              % direction to the observer
g = alpha - phi;
it = find(abs(B) > 1e-3*B0, 1);                     % start searching where the field starts
it = it - 1 + find(g(it:end) >= 0, 1);
R_src = norm(obs - traj.r(:, it)');
fprintf('tangent point: z = %.5f m, field there %.4f T (%.3f of B0), R_src = %.4f m\n', z(it), abs(B(it)), abs(B(it))/B0, R_src);
% the analytic formula needs the arc to extend well beyond the formation angle on both sides
fprintf('formation angle (lambda/rho)^(1/3): %.2f mrad at 15 nm, %.2f mrad at 150 nm; 1/gamma = %.2f mrad\n', ...
    1e3*(15e-9/rho)^(1/3), 1e3*(150e-9/rho)^(1/3), 1e3/gam);

% analytic, theta = 0 (sigma polarisation, in the orbit plane = y here)
yv = w / wc;
S_an = 3*e0^2*gam^2/(16*pi^3*eps0*c0) * yv.^2 .* besselk(2/3, yv/2).^2;
fac = (R_src / norm(obs))^2;
S_mag = Sm_tot * fac;                                  % per sr seen from the tangent point
S_mag_y = Sm_xyz(2, :) * fac;
r = S_mag ./ S_an;
lam_chk = [15 20 30 50 70 100 150];
ichk = arrayfun(@(l) find(abs(lam - l) == min(abs(lam - l)), 1), lam_chk);
fprintf('MAGNETO (on axis) / analytic (theta = 0) at %s nm: %s\n', mat2str(lam_chk), mat2str(r(ichk), 3));
fprintf('min %.3f, max %.3f over 15-150 nm; sigma (y) polarisation share of MAGNETO: %.6f\n', min(r), max(r), min(S_mag_y ./ S_mag));

cmpd.lambda_nm = lam; cmpd.omega = w; cmpd.S_magneto_per_sr_source = S_mag; cmpd.S_analytic = S_an; cmpd.ratio = r;
cmpd.rho_m = rho; cmpd.omega_c = wc; cmpd.B0_T = B0; cmpd.gamma = gam; cmpd.R_src_m = R_src; cmpd.z_tangent_m = z(it);
cmpd.note = 'd2W/domega dOmega [J/(rad/s)/sr] per sr seen from the tangent point; analytic = Jackson (14.79), theta = 0';
save(fullfile(res_dir, 'compare_dipole_onaxis.mat'), 'cmpd');
fid = fopen(fullfile(res_dir, 'compare_dipole_onaxis.txt'), 'wt');
fprintf(fid, '%% 240 MeV dipole, observer on axis at (0,0,11) m; d2W/domega dOmega [J/(rad/s)/sr] per sr seen from the tangent point (z = %.4f m, R_src = %.4f m)\n', z(it), R_src);
fprintf(fid, '%% Columns: 1 lambda [nm]  2 omega [rad/s]  3 MAGNETO  4 analytic Jackson (14.79) theta = 0 (rho = %.4f m, lambda_c = %.3f nm)  5 ratio 3/4\n', rho, 2*pi*c0/wc*1e9);
fprintf(fid, '%-14.8f %-14.6e %-14.6e %-14.6e %-14.6e\n', [lam; w; S_mag; S_an; r]);
fclose(fid);

xt = [15 20 30 50 70 100 150];
f = figure('Position', [100 100 900 650]);
subplot(2,1,1);
loglog(lam, S_mag, 'k-', lam, S_an, 'r--', 'LineWidth', 1.3); grid on; xlim([15 150]); set(gca, 'XTick', xt);
legend('MAGNETO, on axis (fan edge)', 'Jackson (14.79), \theta = 0 (long arc, uniform field)', 'Location', 'southeast');
ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]');
title(sprintf('240 MeV dipole (\\rho = %.3f m, \\lambda_c = %.1f nm): on-axis vs analytic', rho, 2*pi*c0/wc*1e9));
subplot(2,1,2);
semilogx(lam, r, 'k-', 'LineWidth', 1.3); grid on; xlim([15 150]); set(gca, 'XTick', xt);
xlabel('\lambda [nm]'); ylabel('MAGNETO / analytic');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'compare_dipole_onaxis.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'compare_dipole_onaxis.fig'));
