%% ideal_contamination
% Single ideal electron (on axis, no angle, nominal energy; undulator turned over, dz = 60 um): dipole
% contamination of the x-polarised light at the 5 observers, incl. interference, band 50-300 nm.
% The ideal electron leaves the undulator with x' = +4.5e-7 rad (residual steering of the field map, turned over),
% so it is not exactly in the dipole bend plane: this is the physical (not noise) source of dipole x light at x = 0.
% Output: results/ideal_contamination.txt

clear;
this_dir = fileparts(mfilename('fullpath')); root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
c0 = 299792458; eps0 = 8.85418781e-12; lam = logspace(log10(50), log10(300), 300)*1e-9; om = 2*pi*c0./lam;
W = @(S) -trapz(om, S);
tr = read_gpt_trajectory(fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole', 'ideal_electron', 'result_sm_dz60.gdf'));
seg = {sub(tr, tr.r(3,:) < 0.63), sub(tr, tr.r(3,:) >= 0.63)};
fid = fopen(fullfile(this_dir, 'results', 'ideal_contamination.txt'), 'w');
fprintf(fid, '%% single ideal electron, x polarisation, 50-300 nm [J/sr]\n%% y_mm W_und_x W_dip_x W_tot_x contamination=(W_tot-W_und)/W_und dip_over_und\n');
for y = [-1 -1.5 -2 -2.5 -3]*1e-3
    obs = [0 y 5.9]; pref = eps0*c0*norm(obs)^2/pi;
    Eu = emission_time_spectrum(seg{1}, obs, lam, 'fB', norm(obs)/c0).Ew(1,:);
    Ed = emission_time_spectrum(seg{2}, obs, lam, 'fB', norm(obs)/c0).Ew(1,:);
    wu = W(pref*abs(Eu).^2); wd = W(pref*abs(Ed).^2); wt = W(pref*abs(Eu+Ed).^2);
    fprintf(fid, '%5.1f %.4e %.4e %.4e %+.3e %.3e\n', 1e3*y, wu, wd, wt, (wt-wu)/wu, wd/wu);
end
fclose(fid); type(fullfile(this_dir, 'results', 'ideal_contamination.txt'));

function t = sub(t, m)
k = find(m); kf = k(any(t.fB(:, k) ~= 0, 1)); k = max(k(1), kf(1)-1) : min(k(end), kf(end)+1);
for f = {'t', 'G'};             t.(f{1}) = t.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; t.(f{1}) = t.(f{1})(:, k); end
end
