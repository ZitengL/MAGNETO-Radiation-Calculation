%% nokick_check
% Is the x-polarised dipole light of the ideal electron at x = 0 caused by the undulator's residual x' kick?
% Compare the ideal electron (undulator turned over, result_sm_dz60.gdf, leaves the undulator with x' = +4.5e-7)
% with the same electron launched at x' = -4.547e-7 (result_sm_nokick_dz60.gdf, run_nokick.bat), which enters the
% dipole parallel to the bend plane. Band 50-300 nm and three wavelengths, x polarisation.
% Output: results/nokick_check.txt

clear;
this_dir = fileparts(mfilename('fullpath')); root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole', 'ideal_electron');
c0 = 299792458; eps0 = 8.85418781e-12; lam = logspace(log10(50), log10(300), 300)*1e-9; om = 2*pi*c0./lam; lnm = lam*1e9;
W = @(S) -trapz(om, S);
fid = fopen(fullfile(this_dir, 'results', 'nokick_check.txt'), 'w');
for f = {'result_sm_dz60.gdf', 'result_sm_nokick_dz60.gdf'}
    tr = read_gpt_trajectory(fullfile(gdir, f{1}));
    k = find(tr.r(3,:) > 0.625 & tr.r(3,:) < 0.639, 1);
    fprintf(fid, '%% %s: between the magnets x = %.3e m, x'' = %.3e rad\n', f{1}, tr.r(1,k), tr.b(1,k)/tr.b(3,k));
    fprintf(fid, '%% y_mm W_und_x W_dip_x dip/und(band) S_tot/(S_und+S_dip) at 100, 200, 300 nm\n');
    seg = {sub(tr, tr.r(3,:) < 0.63), sub(tr, tr.r(3,:) >= 0.63)};
    for y = [-1 -1.5 -2 -2.5 -3]*1e-3
        obs = [0 y 5.9]; pref = eps0*c0*norm(obs)^2/pi;
        Eu = emission_time_spectrum(seg{1}, obs, lam, 'fB', norm(obs)/c0).Ew(1,:);
        Ed = emission_time_spectrum(seg{2}, obs, lam, 'fB', norm(obs)/c0).Ew(1,:);
        Q = abs(Eu+Ed).^2./(abs(Eu).^2 + abs(Ed).^2);
        [~, w1] = min(abs(lnm-100)); [~, w2] = min(abs(lnm-200)); [~, w3] = min(abs(lnm-300));
        fprintf(fid, '%5.1f %.4e %.4e %.3e %.4f %.4f %.4f\n', 1e3*y, W(pref*abs(Eu).^2), W(pref*abs(Ed).^2), ...
            W(abs(Ed).^2)/W(abs(Eu).^2), Q(w1), Q(w2), Q(w3));
    end
end
fclose(fid); type(fullfile(this_dir, 'results', 'nokick_check.txt'));

function t = sub(t, m)
k = find(m); kf = k(any(t.fB(:, k) ~= 0, 1)); k = max(k(1), kf(1)-1) : min(k(end), kf(end)+1);
for f = {'t', 'G'};             t.(f{1}) = t.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; t.(f{1}) = t.(f{1})(:, k); end
end
