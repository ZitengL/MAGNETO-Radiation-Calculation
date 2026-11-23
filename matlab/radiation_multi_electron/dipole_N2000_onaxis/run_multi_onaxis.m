%% run_multi_onaxis
% MAGNETO multi-particle calculation (no space charge): 2000 electrons through the
% short dipole (240 MeV, 0.1 m, 0.65 T), observer on axis at z = 11 m, 15-150 nm.
% GPT output step 3.75 um (needed for convergence over 15-150 nm; the earlier 60 um run and its
% results are in results_dz60um/, figures_dz60um/).
% Trajectories: gpt_sims/multi_electron/E240MeV_dipole_N2000/result.gdf, read in
% batches of 100 particles (read_gpt_particles.m).
% For each particle the frequency-domain field E~_j(w) is computed with
% emission_time_spectrum.m (fB method) using one COMMON observer-time reference
% t_ref for all particles, so that the fields can be summed coherently:
%   E~_tot(w) = sum_j E~_j(w),    d2W/dw dOmega = eps0 c R^2 |E~_tot|^2 / pi
% Also kept: the incoherent sum  sum_j |E~_j|^2  (same prefactor), and every E~_j
% (single precision) for later analysis.
% Each batch is saved when done; a rerun continues with the remaining batches.
%
% Outputs (relative to this folder):
%   results/multi_onaxis.mat (sums, spectra), results/multi_onaxis_Ej.mat (E~_j, v7.3),
%   results/multi_onaxis.txt

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
% (read_gpt_particles and particle_traj are in magneto/)                                       % read_gpt_particles, particle_traj
res_dir = fullfile(this_dir, 'results'); if ~exist(res_dir, 'dir'); mkdir(res_dir); end

%% 1. Settings
run_dir   = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E240MeV_dipole_N2000');
gdf_file  = fullfile(run_dir, 'result.gdf');
N         = 2000;
batch     = 100;                           % tout 3.75 um: ~75k samples per particle, keep batches small (memory)
obs_point = [0, 0, 11];
lambda_nm = logspace(log10(15), log10(150), 5000);
lambda    = lambda_nm * 1e-9;
c0 = 299792458; eps0 = 8.85418781e-12;
t_ref     = norm(obs_point) / c0;          % common reference: light leaving the origin at t = 0
nb        = ceil(N / batch);
D0 = readmatrix(fullfile(run_dir, 'initial_distribution.txt'), 'NumHeaderLines', 1);   % x y z GBx GBy GBz G nmacro

state_file = fullfile(res_dir, 'multi_onaxis_state.mat');
ej_file    = fullfile(res_dir, 'multi_onaxis_Ej.mat');
if exist(state_file, 'file')
    S = load(state_file); Esum = S.Esum; incoh = S.incoh; done_b = S.done_b;
    fprintf('Continuing: %d of %d batches done\n', nnz(done_b), nb);
else
    Esum = zeros(3, numel(lambda)); incoh = zeros(3, numel(lambda)); done_b = false(1, nb);
    if exist(ej_file, 'file'); delete(ej_file); end
end
mE = matfile(ej_file, 'Writable', true);

%% 2. Batches
t0 = tic;
for b = find(~done_b)
    ids = (b-1)*batch + 1 : min(b*batch, N);
    tb = tic;
    tr = retry_oom(@() read_gpt_particles(gdf_file, ids, true));     % no electric field in the dipole: skip fE
    % particle IDs follow the order of initial_distribution: the transverse positions
    % before the undulator are still the initial ones (drift, no divergence)
    dx = max(abs(tr.x(:, 1) - D0(ids, 1))); dy = max(abs(tr.y(:, 1) - D0(ids, 2)));
    if dx > 1e-12 || dy > 1e-12; error('ID / initial distribution mismatch (%.2e, %.2e)', dx, dy); end
    tread = toc(tb);
    Eb = complex(zeros(numel(ids), 3, numel(lambda), 'single'));
    for k = 1:numel(ids)
        traj = field_region(particle_traj(tr, k));
        s = retry_oom(@() emission_time_spectrum(traj, obs_point, lambda, 'fB', t_ref));
        Esum  = Esum + s.Ew;
        incoh = incoh + abs(s.Ew).^2;
        Eb(k, :, :) = single(s.Ew);
    end
    mE.Ej(ids, 1:3, 1:numel(lambda)) = Eb;
    done_b(b) = true;
    save(state_file, 'Esum', 'incoh', 'done_b');
    fprintf('batch %2d / %d (IDs %d-%d): read %.0f s, total %.0f s; %.1f min elapsed\n', b, nb, ids(1), ids(end), ...
        tread, toc(tb), toc(t0)/60);
    clear tr Eb
end

%% 3. Spectra and save
pref = eps0 * c0 * norm(obs_point)^2 / pi;
comp = {'x', 'y', 'z'};
for c = 1:3
    out.d2W_dwdOmega.(comp{c})      = pref * abs(Esum(c, :)).^2;      % coherent sum, J/(rad/s)/sr
    out.incoherent_sum.(comp{c})    = pref * incoh(c, :);            % sum_j |E_j|^2
end
out.d2W_dwdOmega.total   = out.d2W_dwdOmega.x + out.d2W_dwdOmega.y + out.d2W_dwdOmega.z;
out.incoherent_sum.total = out.incoherent_sum.x + out.incoherent_sum.y + out.incoherent_sum.z;
out.Esum = Esum; out.lambda_nm = lambda_nm; out.omega = 2*pi*c0 ./ lambda;
out.params.gdf_file = gdf_file; out.params.N = N; out.params.obs_point_m = obs_point; out.params.t_ref_s = t_ref;
out.params.note = 'd2W_dwdOmega = coherent sum |sum_j E_j|^2; incoherent_sum = sum_j |E_j|^2; per sr of the sphere R = 11 m centred at the origin';
out.params.date = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
multi = out; %#ok<NASGU>
save(fullfile(res_dir, 'multi_onaxis.mat'), 'multi');
fid = fopen(fullfile(res_dir, 'multi_onaxis.txt'), 'wt');
fprintf(fid, '%% 2000 electrons (no space charge), short dipole (0.1 m, 0.65 T), 240 MeV, on axis at (0, 0, 11) m\n');
fprintf(fid, '%% Columns: 1 lambda [nm]  2 omega [rad/s]  3-6 d2W/domega dOmega of the coherent sum [J/(rad/s)/sr] total,x,y,z  7 incoherent sum sum_j|E_j|^2 (same units), total\n');
fprintf(fid, '%-14.8f %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e %-14.6e\n', [lambda_nm; out.omega; out.d2W_dwdOmega.total; ...
    out.d2W_dwdOmega.x; out.d2W_dwdOmega.y; out.d2W_dwdOmega.z; out.incoherent_sum.total]);
fclose(fid);
fprintf('Finished in %.1f min\n', toc(t0)/60);

function out = retry_oom(fun)
% call fun(); on out-of-memory (other programs can use up the system memory) wait and retry
for attempt = 1:30
    try
        out = fun();
        return
    catch err
        if ~contains(err.message, 'memory', 'IgnoreCase', true) || attempt == 30; rethrow(err); end
        fprintf('out of memory (attempt %d), waiting 60 s\n', attempt);
        pause(60);
    end
end
end

function traj = field_region(traj)
% Keep only samples with nonzero field (plus one zero-field sample on each side).
k = find(any(traj.fB ~= 0, 1));
k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
