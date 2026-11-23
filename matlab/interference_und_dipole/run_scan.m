%% run_scan
% Undulator (real Radia map, centre z = 1 m) + short dipole downstream at distance D, single electron,
% 240 MeV. For every D: GPT run of combined.in and dipole_only.in (gpt_sims/single_electron/
% E240MeV_und15p_dipole_interference, D on the command line), then the frequency-domain field
% E~(w) (x, y, z) with emission_time_spectrum.m (fB) at two observers:
%   'onaxis' (0, 0, 11) m;  'th0p5' 11 m (sin 0.5 mrad, 0, cos 0.5 mrad)  (horizontal = wiggle plane)
% with the common reference time t_ref = |obs|/c, so that fields of different runs can be added.
% The undulator-only field (independent of D) is computed once.
% D list: fine scan 0.500:0.001:0.700 m and coarse scan 0.1:0.1:1.5 m.
% Resumable: one file per D in results/scan/ (skipped if present). GDF files are kept in the GPT folder.
% The dipole-only trajectory is aligned with the dipole part of the combined one (time and x(z)), so that
% E_comb = E_und + E_dip (checked at D = 0.5 m to <= 8.4e-5, check_D0p5.m).
% Output: results/scan/D<mm>mm.mat (Ew_comb, Ew_dip: 2 x 3 x Nw single complex; dt_align;
%         t_arr_entrance = arrival time at the 2 observers of light leaving the dipole entrance with the
%         electron, for the delay prediction), results/undulator_only.mat

clear;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(this_dir));
addpath(fullfile(root_dir, 'magneto'));
gdir = fullfile(root_dir, 'gpt_sims', 'single_electron', 'E240MeV_und15p_dipole_interference');
sdir = fullfile(this_dir, 'results', 'scan'); if ~exist(sdir, 'dir'); mkdir(sdir); end
gpt = fullfile(getenv('GPTBIN'), 'gpt.exe'); lic = getenv('GPTLICENSE');
c0 = 299792458;
lambda_nm = logspace(log10(15), log10(150), 5000); lambda = lambda_nm*1e-9;
obs = {[0 0 11], 11*[sin(0.5e-3) 0 cos(0.5e-3)]}; obs_name = {'onaxis', 'th0p5'};
Dlist = unique(round([0.500:0.001:0.700, 0.1:0.1:1.5]*1e3))/1e3;

%% undulator only (once)
uf = fullfile(this_dir, 'results', 'undulator_only.mat');
if ~exist(uf, 'file')
    T = field_region(read_gpt_trajectory(fullfile(gdir, 'result_undulator_only.gdf')));
    Ew_und = complex(zeros(2, 3, numel(lambda), 'single'));
    for q = 1:2; s = emission_time_spectrum(T, obs{q}, lambda, 'fB', norm(obs{q})/c0); Ew_und(q, :, :) = single(s.Ew); end
    save(uf, 'Ew_und', 'lambda_nm', 'obs', 'obs_name');
end

%% scan
t0 = tic; old = cd(gdir); cleaner = onCleanup(@() cd(old));
for i = 1:numel(Dlist)
    D = Dlist(i); tag = sprintf('D%04dmm', round(D*1e3)); f = fullfile(sdir, [tag '.mat']);
    if exist(f, 'file'); continue; end
    gc = fullfile(gdir, 'runs_combined', ['result_combined_' tag '.gdf']);
    gd = fullfile(gdir, 'runs_dipole_only', ['result_dipole_only_' tag '.gdf']);
    st = system(sprintf('"%s" -o "%s" combined.in D=%.6f GPTLICENSE=%s', gpt, gc, D, lic));
    st = st + system(sprintf('"%s" -o "%s" dipole_only.in D=%.6f GPTLICENSE=%s', gpt, gd, D, lic));
    if st ~= 0; error('GPT failed for D = %.3f', D); end
    Tc = field_region(read_gpt_trajectory(gc)); Td = field_region(read_gpt_trajectory(gd));
    % Align the dipole-only electron with the combined one: in the combined run it reaches the dipole
    % later (slippage in the undulator, c dt ~ 0.5 um) and with the small transverse offset/angle left by
    % the real undulator (~1 um, ~1 urad). Time shift at the dipole centre; x(z), beta_x(z) of the combined run.
    zc = 1.32004 + D + 0.05;
    dt_align = interp1(Tc.r(3,:), Tc.t, zc) - interp1(Td.r(3,:), Td.t, zc);
    Td.t = Td.t + dt_align;
    Td.r(1,:) = interp1(Tc.r(3,:), Tc.r(1,:), Td.r(3,:), 'linear', 'extrap');
    Td.b(1,:) = interp1(Tc.r(3,:), Tc.b(1,:), Td.r(3,:), 'linear', 'extrap');
    % arrival time at the observers of light leaving the dipole entrance (z = zc - 0.05) with the electron
    ze = zc - 0.05; te = interp1(Tc.r(3,:), Tc.t, ze); re = [interp1(Tc.r(3,:), Tc.r(1,:), ze); interp1(Tc.r(3,:), Tc.r(2,:), ze); ze];
    t_arr_entrance = [te + norm(obs{1}(:) - re)/c0, te + norm(obs{2}(:) - re)/c0];
    Ew_comb = complex(zeros(2, 3, numel(lambda), 'single')); Ew_dip = Ew_comb;
    for q = 1:2
        s = emission_time_spectrum(Tc, obs{q}, lambda, 'fB', norm(obs{q})/c0); Ew_comb(q, :, :) = single(s.Ew);
        s = emission_time_spectrum(Td, obs{q}, lambda, 'fB', norm(obs{q})/c0); Ew_dip(q, :, :)  = single(s.Ew);
    end
    exit_angle = atan2(Tc.b(2, end), Tc.b(3, end));
    save(f, 'D', 'Ew_comb', 'Ew_dip', 'exit_angle', 'dt_align', 't_arr_entrance');
    fprintf('%3d / %d: D = %.3f m done, %.1f min elapsed\n', i, numel(Dlist), D, toc(t0)/60);
end
fprintf('Scan finished in %.1f min\n', toc(t0)/60);

function traj = field_region(traj)
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
