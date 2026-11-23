%% diag_textbook
% Diagnostic (not a result for the paper): is the difference between MAGNETO on the unduplan
% trajectory and the analytic (26.46) a property of the method, or of the different end of the field?
% The analytic formula assumes the textbook trajectory: exactly N = 15 periods, beta_x = (K/gamma) sin(ku z'),
% beta_x = 0 at both ends, straight outside. This is the exact motion in a HARD-EDGED cosine field
%   By(z) = B0 cos(ku z'),  z' = z - 1 m,  |z'| <= N lambda_u/2,  0 outside   (B0 = K 2 pi m c/(e lambda_u))
% (it starts and ends at -B0; unduplan instead has 1/4, 3/4 end half periods and a sine phase).
% Exact motion (no tracking): p_x = e int By dz -> gamma beta_x = K sin(ku z'), gamma constant,
% beta_z = sqrt(beta^2 - beta_x^2), t(z) = int dz/(c beta_z), x(z) = int beta_x/beta_z dz (cumtrapz, 0.25 um).
% Samples every dzs in z, with sample points exactly at both field edges. MAGNETO = emission_time_spectrum.m (fB).
% Compared with (26.46) at the same 4 observers as compare_undulator.m.
clear; close all;
this_dir = fileparts(mfilename('fullpath'));
up_dir = fileparts(this_dir); root_dir = fileparts(fileparts(fileparts(up_dir)));
addpath(fullfile(root_dir, 'magneto')); addpath(up_dir);
c0 = 299792458; me_g = 9.10953e-31; qe_g = 1.6021892e-19;    % GPT constants, as in the GPT runs
gam = 470.66442480; K = 0.980309; lamu = 0.030472767; Np = 15; zu = 1.0; ku = 2*pi/lamu;
beta = sqrt(1 - 1/gam^2); B0 = K*2*pi*me_g*c0/(qe_g*lamu);
lam_nm = logspace(log10(15), log10(150), 5000); lam = lam_nm*1e-9;
nw = round(0.01 / mean(diff(log(lam_nm)))); env = @(S) movmax(S, 2*nw + 1);

Lh = Np*lamu/2;
zf = linspace(-Lh, Lh, round(2*Lh/0.25e-6) + 1);              % fine grid for the quadrature
bx = K/gam * sin(ku*zf); bz = sqrt(beta^2 - bx.^2);
tf = cumtrapz(zf, 1 ./ (c0*bz)); xf = cumtrapz(zf, bx ./ bz);
thu5 = atan(11*sin(0.5e-3) / (11*cos(0.5e-3) - 1));
obs = {[0 0 11], 11*[sin(0.5e-3) 0 cos(0.5e-3)], [0 0 zu+1000], [0 0 zu] + 1000*[sin(thu5) 0 cos(thu5)]};
thu = [0 thu5 0 thu5]; names = {'11m_theta0', '11m_theta0p5', '1km_theta0', '1km_theta0p5'};

for dzs = [7.5e-6 3.75e-6]
    step = round(dzs/0.25e-6); k = 1:step:numel(zf); if k(end) ~= numel(zf); k = [k numel(zf)]; end
    tr.t = tf(k); tr.r = [xf(k); zeros(1, numel(k)); zu + zf(k)]; tr.b = [bx(k); zeros(1, numel(k)); bz(k)];
    tr.G = gam*ones(1, numel(k)); tr.fE = zeros(3, numel(k)); tr.fB = [zeros(1, numel(k)); B0*cos(ku*zf(k)); zeros(1, numel(k))];   % electron: dp_x/dz = e By
    tr.q = -qe_g; tr.m = me_g; tr.nmacro = 1;
    % sign check: with this By the Lorentz force must give d(beta_x)/dt of the prescribed beta_x
    bd = (tr.q/(gam*tr.m*c0)) * c0 * (tr.b(3,:).*tr.fB(2,:) * -1);   % x component of q/(g m c) c (b x B), b x B _x = by Bz - bz By
    bd_num = gradient(tr.b(1,:), tr.t);
    fprintf('\nsample step %.2f um (%d samples); max |dbx/dt (Lorentz) - dbx/dt (numerical)| / max = %.1e\n', dzs*1e6, numel(k), max(abs(bd - bd_num))/max(abs(bd_num)));
    for q = 1:4
        a = undulator_analytic(lam, gam, K, Np, lamu, thu(q), 0);
        s = emission_time_spectrum(tr, obs{q}, lam, 'fB');
        St = s.d2W_dwdOmega.total * (norm(obs{q} - [0 0 zu]) / norm(obs{q}))^2;
        e = abs(St - a.W46.total) ./ env(a.W46.total);
        sel = a.W46.total > 1e-2*max(a.W46.total);
        fprintf('  %-13s |textbook-trajectory MAGNETO - (26.46)| / envelope: max %.2e, median %.1e; where S > 1%% of max: ratio - 1 in [%+.2e, %+.2e]\n', ...
            names{q}, max(e), median(e), min(St(sel)./a.W46.total(sel) - 1), max(St(sel)./a.W46.total(sel) - 1));
        D.(['dz' num2str(round(dzs*1e8))]).(['c_' names{q}]) = struct('S_textbook_traj', St, 'S_analytic', a.W46.total, 'err_env', e);
        if dzs < 5e-6
            if thu(q) == 0; kl = [1 3 5]; else; kl = 1:5; end
            lk = a.lambda1*1e9 ./ kl; r = zeros(size(kl));
            for j = 1:numel(kl)
                win = abs(lam_nm - lk(j)) < 0.04*lk(j);
                [~, va] = peak_refined(lam_nm, a.W46.total, win); [~, vt] = peak_refined(lam_nm, St, win); r(j) = vt/va;
            end
            fprintf('      peak value textbook-trajectory MAGNETO / (26.46), k = %s: %s\n', mat2str(kl), mat2str(r, 4));
            f = figure('Position', [100 100 1000 650]);
            subplot(3,1,[1 2]); loglog(lam_nm, a.W46.total, 'k-', lam_nm, St, 'g--', 'LineWidth', 1.2); grid on; xlim([15 150]);
            set(gca, 'XTick', [15 20 30 40 50 70 100 150]); ylim(max(a.W46.total)*[1e-6 3]);
            legend('analytic (26.46), N = 15', 'MAGNETO on the textbook trajectory (hard-edged cos field)', 'Location', 'southeast');
            ylabel('d^2W/d\omega d\Omega [J/(rad/s)/sr]'); title(['Diagnostic: textbook trajectory, ' strrep(names{q}, '_', ' ')]);
            subplot(3,1,3); semilogy(lam_nm, e, 'g-'); set(gca, 'XScale', 'log'); grid on; xlim([15 150]); ylim([1e-6 10]);
            set(gca, 'XTick', [15 20 30 40 50 70 100 150]); xlabel('\lambda [nm]'); ylabel('diff / local envelope');
            set(findall(f, 'Type', 'axes'), 'Toolbar', []);
            if ~exist(fullfile(this_dir, 'figures'), 'dir'); mkdir(fullfile(this_dir, 'figures')); end
            exportgraphics(f, fullfile(this_dir, 'figures', ['diag_textbook_' names{q} '.png']), 'Resolution', 150); close(f);
        end
    end
end
D.lambda_nm = lam_nm; D.note = 'MAGNETO on the exact textbook trajectory (hard-edged cosine field, N = 15) vs (26.46); per sr from the undulator centre';
if ~exist(fullfile(this_dir, 'results'), 'dir'); mkdir(fullfile(this_dir, 'results')); end
save(fullfile(this_dir, 'results', 'diag_textbook.mat'), 'D');

function [lp, Sp] = peak_refined(lam, S, win)
[~, j] = max(S .* win);
x = log(lam(j-1:j+1)); y = log(S(j-1:j+1)); p = polyfit(x - x(2), y, 2);
dx = -p(2)/(2*p(1)); lp = exp(x(2) + dx); Sp = exp(polyval(p, dx));
end
