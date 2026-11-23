%% check_unduplan
% Checks of the ideal-undulator run (unduplan, Nu = 16, trims 1/4, 3/4) and comparison of its
% harmonic wavelengths with the real (Radia) 15-period map run (E240MeV_und15p_paraxial):
%   1. field: extent, peak, taper of the end half periods, deviation from Bu*Fs*sin(ku*(z - zund))
%   2. trajectory: exit angle and offset, gamma
%   3. resonant wavelengths: <beta_z> method (harmonic_wavelengths.m) and formula
%      lambda_n = lambda_u/(2 n gamma^2)(1 + K^2/2) with the input K, for both runs
%   4. peak positions of the computed spectra near each harmonic (on axis n = 1, 3, 5;
%      theta_x = 0.5 mrad on the 11 m sphere n = 1..5), ideal vs real map, 15-150 nm, 5000 log points,
%      emission_time_spectrum.m (fB); tout convergence of the ideal run (60 vs 15 um)
clear;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
c0 = 299792458;
zund = 1.0; Nu = 16; lamu = 0.030472767; Kpk = 0.980309; ku = 2*pi/lamu;
me_g = 9.10953e-31; qe_g = 1.6021892e-19;                 % GPT built-in constants (used in the .in)
Bu = Kpk*2*pi*me_g*c0/(qe_g*lamu);

id = read_gpt_trajectory(fullfile(this_dir, 'result.gdf'));
z = id.r(3, :); k = find(any(id.fB ~= 0, 1));
fprintf('Bu = %.6f T; field nonzero for z = %.5f .. %.5f m (expected %.5f .. %.5f)\n', Bu, z(k(1)), z(k(end)), zund - Nu*lamu/2, zund + Nu*lamu/2);
fprintf('max |fBx| %.2e T, max |fBy| %.6f T, max |fBz| %.2e T (on axis: Bz = sinh(ku y) part, y = 0)\n', max(abs(id.fB(1,:))), max(abs(id.fB(2,:))), max(abs(id.fB(3,:))));
% half-period peaks
z0 = zund - Nu*lamu/2; hp = lamu/2; pk = zeros(1, 2*Nu);
for m = 1:2*Nu
    sel = z >= z0 + (m-1)*hp & z < z0 + m*hp; pk(m) = max(abs(id.fB(2, sel)));
end
fprintf('half-period peaks / Bu: first 3 %s, last 3 %s, middle min/max %.6f / %.6f\n', mat2str(pk(1:3)/Bu, 4), mat2str(pk(end-2:end)/Bu, 4), min(pk(3:end-2))/Bu, max(pk(3:end-2))/Bu);
Fs = ones(size(z)); m = floor((z - z0)/hp) + 1;
Fs(m == 1 | m == 2*Nu) = 0.25; Fs(m == 2 | m == 2*Nu-1) = 0.75; Fs(m < 1 | m > 2*Nu) = 0;
Bid = Bu * Fs .* sin(ku*(z - zund));
fprintf('max |fBy - Bu Fs sin(ku (z - zund))| = %.2e T (sign of the sine: %+d)\n', min(max(abs(id.fB(2,:) - Bid)), max(abs(id.fB(2,:) + Bid))), sign(-max(abs(id.fB(2,:) + Bid)) + max(abs(id.fB(2,:) - Bid))) * -1);
fprintf('gamma %.8f .. %.8f; exit beta_x %.3e, exit x %.3e m, mean x over the field %.3e m, max |x| %.3f um, K/gamma %.4e\n', ...
    min(id.G), max(id.G), id.b(1,end), id.r(1,end), mean(id.r(1,k)), 1e6*max(abs(id.r(1,:))), Kpk/id.G(1));

%% resonant wavelengths
rm = read_gpt_trajectory(fullfile(root_dir, 'gpt_sims', 'single_electron', 'E240MeV_und15p_paraxial', 'result.gdf'));
hi = harmonic_wavelengths(id, 1:5); hr = harmonic_wavelengths(rm, 1:5);
gam = id.G(1); lam_f = lamu/(2*gam^2)*(1 + Kpk^2/2) ./ (1:5);
fprintf('\n<beta_z> method:  ideal K %.6f, lambda_u %.6f mm | real map K %.6f, lambda_u %.6f mm\n', hi.K_from_beta, hi.lambda_u*1e3, hr.K_from_beta, hr.lambda_u*1e3);
fprintf('lambda_n [nm]    n=1        n=3        n=5\n');
fprintf('formula (K in) %9.4f %10.4f %10.4f\n', lam_f([1 3 5])*1e9);
fprintf('ideal  <beta_z> %9.4f %10.4f %10.4f\n', hi.lambda_beta([1 3 5])*1e9);
fprintf('real   <beta_z> %9.4f %10.4f %10.4f\n', hr.lambda_beta([1 3 5])*1e9);

%% spectra: peak positions
fr = @(t) t;                                         % field region
lam = logspace(log10(15), log10(150), 5000);
idf = field_region(id); rmf = field_region(rm);
i15 = field_region(read_gpt_trajectory(fullfile(this_dir, 'tout_convergence', 'result_dz15um.gdf')));
R = 11; th = [0 0.5e-3]; nw = round(0.01 / mean(diff(log(lam)))); env = @(S) movmax(S, 2*nw + 1);
thu = atan(R*sin(th) ./ (R*cos(th) - zund));
for a = 1:2
    D = R*[sin(th(a)) 0 cos(th(a))];
    Si = emission_time_spectrum(idf, D, lam*1e-9, 'fB'); Si = Si.d2W_dwdOmega.total;
    Sr = emission_time_spectrum(rmf, D, lam*1e-9, 'fB'); Sr = Sr.d2W_dwdOmega.total;
    S15 = emission_time_spectrum(i15, D, lam*1e-9, 'fB'); S15 = S15.d2W_dwdOmega.total;
    fprintf('\ntheta_x = %.1f mrad (theta_u = %.4f mrad): tout 60 vs 15 um max |diff|/env = %.2e\n', th(a)*1e3, thu(a)*1e3, max(abs(Si - S15)./env(S15)));
    if a == 1; nl = [1 3 5]; else; nl = 1:5; end
    lam_n = lamu/(2*gam^2)*(1 + Kpk^2/2 + gam^2*thu(a)^2) ./ nl * 1e9;
    fprintf('  n   formula [nm]   ideal peak [nm]   real-map peak [nm]   ideal/real - 1   peak S ideal/real   ideal peak / max S at this angle\n');
    for q = 1:numel(nl)
        win = abs(lam - lam_n(q)) < 0.04*lam_n(q);
        [li, mi] = peak_refined(lam, Si, win); [lr, mr] = peak_refined(lam, Sr, win);
        fprintf('  %d  %10.4f   %12.4f   %14.4f   %+12.2e   %10.4f   %12.2e\n', nl(q), lam_n(q), li, lr, li/lr - 1, mi/mr, mi/max(Si));
    end
    f = figure('Position', [100 100 950 500]);
    loglog(lam, Si, 'b-', lam, Sr, 'r--', 'LineWidth', 1.1); hold on; grid on; xlim([15 150]);
    for q = 1:numel(nl); xline(lam_n(q), 'k:', sprintf('n = %d', nl(q))); end
    set(gca, 'XTick', [15 20 30 40 50 70 100 150]);
    legend('ideal (unduplan, Nu = 16)', 'real map (Radia, paraxial)', 'Location', 'southeast');
    xlabel('\lambda [nm]'); ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]');
    title(sprintf('240 MeV single electron, \\theta_x = %.1f mrad on the 11 m sphere: ideal vs real undulator', th(a)*1e3));
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    if ~exist(fullfile(this_dir, 'figures'), 'dir'); mkdir(fullfile(this_dir, 'figures')); end
    exportgraphics(f, fullfile(this_dir, 'figures', sprintf('check_unduplan_theta%d.png', a-1)), 'Resolution', 150); close(f);
end

function [lp, Sp] = peak_refined(lam, S, win)
% maximum inside win, refined by a parabola through the 3 points around it (in log lambda, log S)
[~, j] = max(S .* win);
x = log(lam(j-1:j+1)); y = log(S(j-1:j+1)); p = polyfit(x - x(2), y, 2);
dx = -p(2)/(2*p(1)); lp = exp(x(2) + dx); Sp = exp(polyval(p, dx));
end

function traj = field_region(traj)
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
