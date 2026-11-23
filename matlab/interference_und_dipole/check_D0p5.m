%% check_D0p5
% Checks at D = 0.5 m before the scan (runs from gpt_sims/.../E240MeV_und15p_dipole_interference/run.bat):
%   1. tout convergence: combined run with output steps 60/3.75 um vs 30/1.875 um (tout_convergence/)
%   2. decomposition: E_combined vs E_undulator_only + E_dipole_only (same t_ref)
%   3. polarisation content of the two sources and of the sum; 45-degree projections, S2, S3
%   4. delay: tau = d arg(Ex Ey*)/d omega (fit), compared with the arrival-time difference of the
%      two sources (from the trajectory, see below)
clear;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(this_dir));
addpath(fullfile(root_dir, 'magneto'));
gdir = fullfile(root_dir, 'gpt_sims', 'single_electron', 'E240MeV_und15p_dipole_interference');
c0 = 299792458; eps0 = 8.85418781e-12;
lam_nm = logspace(log10(15), log10(150), 5000); lam = lam_nm*1e-9; w = 2*pi*c0./lam;
obs = {[0 0 11], 11*[sin(0.5e-3) 0 cos(0.5e-3)]}; on = {'on axis', '0.5 mrad'};
nw = round(0.01 / mean(diff(log(lam_nm)))); env = @(S) movmax(S, 2*nw + 1);
T.c = field_region(read_gpt_trajectory(fullfile(gdir, 'result_combined_D0p5.gdf')));
T.f = field_region(read_gpt_trajectory(fullfile(gdir, 'tout_convergence', 'result_combined_fine_D0p5.gdf')));
T.d = field_region(read_gpt_trajectory(fullfile(gdir, 'result_dipole_only_D0p5.gdf')));
T.u = field_region(read_gpt_trajectory(fullfile(gdir, 'result_undulator_only.gdf')));
% the electron of the combined run reaches the dipole later (slippage in the undulator): align the
% dipole-only trajectory in time (and position) at the dipole centre z = 1.87 m
zc = 1.32004 + 0.5 + 0.05;
dt_al = interp1(T.c.r(3,:), T.c.t, zc) - interp1(T.d.r(3,:), T.d.t, zc);
dx_al = interp1(T.c.r(3,:), T.c.r(1,:), zc);
fprintf('alignment of the dipole-only run: delay %.4e s (c dt = %.4f um), x offset of the combined run at the dipole %.2e m\n', dt_al, dt_al*c0*1e6, dx_al);
T.d.t = T.d.t + dt_al;
T.d.r(1,:) = interp1(T.c.r(3,:), T.c.r(1,:), T.d.r(3,:), 'linear', 'extrap');   % transverse drift of the combined run
T.d.b(1,:) = interp1(T.c.r(3,:), T.c.b(1,:), T.d.r(3,:), 'linear', 'extrap');
for q = 1:2
    O = obs{q}; pref = eps0*c0*norm(O)^2/pi; tr = norm(O)/c0;
    E = struct(); for k = 'cfdu'; s = emission_time_spectrum(T.(k), O, lam, 'fB', tr); E.(k) = s.Ew; end
    S = @(e) pref*sum(abs(e).^2, 1);
    fprintf('\n=== %s\n', on{q});
    fprintf('1. tout: |S(60/3.75) - S(30/1.875)|/env max %.1e; S2: |diff|/env(S0) max %.1e\n', ...
        max(abs(S(E.c) - S(E.f))./env(S(E.f))), max(abs(2*real(E.c(1,:).*conj(E.c(2,:))) - 2*real(E.f(1,:).*conj(E.f(2,:))))*pref ./ env(S(E.f))));
    fprintf('2. |E_comb - (E_und + E_dip)| / |E_comb| (per component, envelope): x %.1e, y %.1e\n', ...
        max(abs(E.c(1,:) - E.u(1,:) - E.d(1,:)) ./ movmax(abs(E.c(1,:)), 2*nw+1)), max(abs(E.c(2,:) - E.u(2,:) - E.d(2,:)) ./ movmax(abs(E.c(2,:)), 2*nw+1)));
    Su = S(E.u); Sd = S(E.d); Sc = S(E.c);
    fprintf('3. undulator: x/total %.4f..%.4f (y/total max %.1e); dipole: y/total %.4f..%.4f (x/total max %.1e)\n', ...
        min(pref*abs(E.u(1,:)).^2./Su), max(pref*abs(E.u(1,:)).^2./Su), max(pref*abs(E.u(2,:)).^2./Su), ...
        min(pref*abs(E.d(2,:)).^2./Sd), max(pref*abs(E.d(2,:)).^2./Sd), max(pref*abs(E.d(1,:)).^2./Sd));
    fprintf('   S_dip/S_und: min %.1e, median %.1e, max %.1e; |S_comb - S_und - S_dip|/S_comb max %.1e (total-intensity interference)\n', ...
        min(Sd./Su), median(Sd./Su), max(Sd./Su), max(abs(Sc - Su - Sd)./Sc));
    S2 = 2*pref*real(E.c(1,:).*conj(E.c(2,:))); S3 = -2*pref*imag(E.c(1,:).*conj(E.c(2,:)));
    V = sqrt(S2.^2 + S3.^2) ./ Sc;                       % = fringe visibility of the 45-degree projection
    fprintf('   45-degree projection: fringe visibility sqrt(S2^2+S3^2)/S0 median %.2e, max %.2e\n', median(V), max(V));
    % 4. delay from the phase of Ex Ey*
    ph = unwrap(angle(E.c(1,:).*conj(E.c(2,:))));
    % only where both components are strong (avoids the pi jumps at the undulator spectral nulls)
    m = abs(E.c(1,:)) > 0.1*movmax(abs(E.c(1,:)), 2*nw+1) & abs(E.c(2,:)) > 0.1*movmax(abs(E.c(2,:)), 2*nw+1);
    fprintf('4. fraction of the grid where both |Ex| and |Ey| > 10%% of their envelope: %.2f\n', mean(m));
end

function traj = field_region(traj)
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
