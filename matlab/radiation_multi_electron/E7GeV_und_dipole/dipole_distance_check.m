% compare GPT dipole at +0.3 m with the translated trajectory of dipole_distance_scan.m
root_dir = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(fullfile(root_dir, 'magneto'));
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole', 'ideal_electron');
c0 = 299792458; lam = logspace(log10(50), log10(300), 300)*1e-9; obs = [0 -2e-3 5.9];
sub = @(t, m) subf(t, m);
t0 = read_gpt_trajectory(fullfile(gdir, 'result_sm_dz60.gdf')); t3 = read_gpt_trajectory(fullfile(gdir, 'result_sm_dL0p3_dz60.gdf'));
u0 = subf(t0, t0.r(3,:) < 0.63); d0 = subf(t0, t0.r(3,:) >= 0.63); d3 = subf(t3, t3.r(3,:) >= 0.63);
v0 = c0*d0.b(:,1); dt = 0.3/v0(3); s = d0; s.t = s.t + dt; s.r = s.r + v0*dt;
Eu = emission_time_spectrum(u0, obs, lam, 'fB', norm(obs)/c0).Ew;
Et = emission_time_spectrum(s, obs, lam, 'fB', norm(obs)/c0).Ew;
Eg = emission_time_spectrum(d3, obs, lam, 'fB', norm(obs)/c0).Ew;
St = sum(abs(Eu + Et).^2, 1); Sg = sum(abs(Eu + Eg).^2, 1);
fprintf('dipole part: max|E_transl - E_GPT| / max|E_GPT| = %.2e\n', max(abs(Et - Eg), [], 'all')/max(abs(Eg), [], 'all'));
fprintf('total spectrum: max |S_transl - S_GPT| / max S = %.2e\n', max(abs(St - Sg))/max(Sg));
function t = subf(t, m)
k = find(m); kf = k(any(t.fB(:, k) ~= 0, 1)); k = max(k(1), kf(1)-1) : min(k(end), kf(end)+1);
for f = {'t', 'G'};             t.(f{1}) = t.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; t.(f{1}) = t.(f{1})(:, k); end
end
