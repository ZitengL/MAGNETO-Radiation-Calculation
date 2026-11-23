%% single_electron
% Ideal electron through the SLAC gap-6 undulator (centre 0.4 m) + dipole (centre 1.3144 m, L = 1 m, 0.13 T, bent up,
% Enge gap 35 mm), GPT runs in gpt_sims/multi_electron/E7GeV_und_dipole/ideal_electron/.
%   sp = undulator as before, sm = undulator turned over (180 deg about z).
% Observers: b_new = (0, -1 mm, 5.9 m) (mirror of point (b) about the zx plane, side opposite to the bend),
%            b_old = (0, +1 mm, 5.9 m) (bend side, for comparison).
% The emission-time integral is split into the undulator part (z < 0.63 m) and the dipole part (z >= 0.63 m);
% between 0.6238 and 0.6394 m the field is exactly 0, so the split is exact and E_total = E_und + E_dip.
% Output: results/single_electron.mat, results/single_electron.txt, figures/single_electron.png

clear; close all;
this_dir = fileparts(mfilename('fullpath')); root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole', 'ideal_electron');
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
if ~exist(res_dir, 'dir'); mkdir(res_dir); end; if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
c0 = 299792458; eps0 = 8.85418781e-12;
lam = logspace(log10(50), log10(300), 300)*1e-9; om = 2*pi*c0./lam;
obsl = {[0 -1e-3 5.9], [0 1e-3 5.9]}; obsn = {'b_new (0,-1mm)', 'b_old (0,+1mm)'};
runs = {'sp_dz60', 'sm_dz60', 'sm_dz30'}; zsplit = 0.63;

for ir = 1:numel(runs)
    tr = read_gpt_trajectory(fullfile(gdir, ['result_' runs{ir} '.gdf']));
    if ir == 1
        % geometry check
        z = tr.r(3,:); k2 = find(z > 2.0, 1);
        yp = tr.b(2,end)/tr.b(3,end); xp = tr.b(1,end)/tr.b(3,end);
        y59 = tr.r(2,end) + yp*(5.9 - z(end)); x59 = tr.r(1,end) + xp*(5.9 - z(end));
        iu = find(z > 0.63, 1); yu = tr.b(2,iu)/tr.b(3,iu); xu = tr.b(1,iu)/tr.b(3,iu);
        fprintf('after undulator (z = %.3f m): x %.3e m, y %.3e m, x'' %.3e, y'' %.3e rad\n', z(iu), tr.r(1,iu), tr.r(2,iu), xu, yu);
        fprintf('after dipole: y'' = %.6f mrad, x'' = %.2e rad; extrapolated to z = 5.9 m: y = %.3f mm, x = %.3e m\n', 1e3*yp, xp, 1e3*y59, x59);
        nz = find(any(tr.fB ~= 0, 1)); zf = z(nz); gapz = zf(find(diff(zf) > 1e-3)) ; % fB = 0 between the two magnets
        fprintf('field-free gap between the magnets: from z = %.4f m to %.4f m\n', gapz(1), zf(find(zf > gapz(1) + 1e-3, 1)));
        fprintf('|fB| at z = 0.63 m: %.2e T\n', norm(tr.fB(:, find(z > 0.63, 1))));
    end
    for io = 1:2
        obs = obsl{io}; tref = norm(obs)/c0;
        seg = {tr.r(3,:) < zsplit, tr.r(3,:) >= zsplit};
        for is = 1:2
            s = emission_time_spectrum(sub(tr, seg{is}), obs, lam, 'fB', tref);
            E{ir, io, is} = s.Ew;                                                    %#ok<SAGROW>
        end
    end
end
pref = @(obs) eps0*c0*norm(obs)^2/pi;
W = @(S) -trapz(om, S);
fprintf('\nsingle electron, energy in 50-300 nm [J/sr] (x / y polarisation / total)\n');
T = [];
for ir = 1:numel(runs)
    for io = 1:2
        p = pref(obsl{io}); Eu = E{ir,io,1}; Ed = E{ir,io,2}; Et = Eu + Ed;
        Su = p*abs(Eu).^2; Sd = p*abs(Ed).^2; St = p*abs(Et).^2;
        fprintf('%s %-15s und %.3e / %.3e / %.3e | dip %.3e / %.3e / %.3e | total %.3e / %.3e / %.3e | dip/und (total pol.) %.3e\n', ...
            runs{ir}, obsn{io}, W(Su(1,:)), W(Su(2,:)), W(sum(Su)), W(Sd(1,:)), W(Sd(2,:)), W(sum(Sd)), W(St(1,:)), W(St(2,:)), W(sum(St)), W(sum(Sd))/W(sum(Su)));
        R.(runs{ir}){io} = struct('obs', obsl{io}, 'Eu', Eu, 'Ed', Ed);
    end
end
% tout check (sm, both observers, each part)
for io = 1:2
    for is = 1:2
        A = abs(E{2,io,is}).^2; B = abs(E{3,io,is}).^2;
        segn = {'und', 'dip'};
        fprintf('tout 60 vs 30 um, %s, %s: max |dS|/max S %.2e, band energy %.2e\n', obsn{io}, segn{is}, ...
            max(abs(sum(A) - sum(B)))/max(sum(B)), trapz(sum(A))/trapz(sum(B)) - 1);
    end
end
% symmetry check: turned-over undulator at b_new = 180 deg rotation of the original undulator at b_old: (Ex,Ey,Ez) -> (-Ex,-Ey,Ez)
d = max(abs(E{2,1,1} - diag([-1 -1 1])*E{1,2,1}), [], 'all')/max(abs(E{1,2,1}), [], 'all');
fprintf('undulator part: sm at b_new vs rotated sp at b_old: max diff / max %.2e\n', d);
d = max(abs(E{2,1,1}(2,:) - E{1,1,1}(2,:)))/max(abs(E{1,1,1}(2,:)));
fprintf('undulator part at b_new, Ey: sm vs sp max diff / max %.2e (Ex changes sign only)\n', d);
% interference term in the 50-300 nm energy, b_new, sm
p = pref(obsl{1}); Eu = E{2,1,1}; Ed = E{2,1,2};
fprintf('b_new, sm: W_und %.3e, W_dip %.3e, interference 2Re(Eu*Ed) %.3e, total %.3e J/sr\n', W(p*sum(abs(Eu).^2)), W(p*sum(abs(Ed).^2)), ...
    W(2*p*real(sum(conj(Eu).*Ed))), W(p*sum(abs(Eu + Ed).^2)));

save(fullfile(res_dir, 'single_electron.mat'), 'R', 'lam', 'om', 'obsl', 'runs');
% txt: sm_dz60 at b_new
p = pref(obsl{1}); Eu = E{2,1,1}; Ed = E{2,1,2}; Et = Eu + Ed;
M = [lam(:)*1e9, om(:), (p*abs(Eu).^2)', (p*abs(Ed).^2)', (p*abs(Et).^2)'];
fid = fopen(fullfile(res_dir, 'single_electron.txt'), 'wt');
fprintf(fid, '%% ideal electron, undulator turned over + dipole bent up, observer (0, -1 mm, 5.9 m), d2W/dwdOmega [J/(rad/s)/sr]\n');
fprintf(fid, '%% lambda_nm omega und_x und_y und_z dip_x dip_y dip_z total_x total_y total_z\n');
fprintf(fid, [repmat('%.10e ', 1, 11) '\n'], M'); fclose(fid);

f = figure('Position', [50 50 1100 420]);
for io = 1:2
    subplot(1,2,io); hold on; p = pref(obsl{io}); Eu = E{2,io,1}; Ed = E{2,io,2}; Et = Eu + Ed;
    plot(lam*1e9, p*sum(abs(Eu).^2), 'b-', 'DisplayName', 'undulator'); plot(lam*1e9, p*sum(abs(Ed).^2), 'r-', 'DisplayName', 'dipole');
    plot(lam*1e9, p*sum(abs(Et).^2), 'k-', 'DisplayName', 'undulator + dipole');
    plot(lam*1e9, p*abs(Ed(1,:)).^2, 'r:', 'DisplayName', 'dipole, x pol.'); plot(lam*1e9, p*abs(Et(1,:)).^2, 'k:', 'DisplayName', 'total, x pol.');
    set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; xlim([50 300]); xlabel('\lambda [nm]'); ylabel('d^2W/d\omegad\Omega [J/(rad/s)/sr]');
    title(sprintf('single electron, %s', obsn{io})); legend('Location', 'best');
end
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'single_electron.png'), 'Resolution', 200); close(f);

function t = sub(t, m)
k = find(m); kf = k(any(t.fB(:, k) ~= 0, 1)); k = max(k(1), kf(1)-1) : min(k(end), kf(end)+1);
for f = {'t', 'G'};             t.(f{1}) = t.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; t.(f{1}) = t.(f{1})(:, k); end
end
