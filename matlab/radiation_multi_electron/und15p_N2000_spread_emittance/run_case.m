function run_case(case_dir, gdf_name, out_file)
%% run_case
% MAGNETO multi-particle calculation (no space charge) for one case of the energy-spread / emittance
% study: 2000 electrons, 15-period undulator (Radia map, +-2 mm Laplace map), 240 MeV, observers
%   'onaxis' (0, 0, 11) m and 'th0p5' 11 m (sin 0.5 mrad, 0, cos 0.5 mrad)  (horizontal = wiggle plane),
% 15-150 nm, 5000 log points, common reference time t_ref = |obs|/c (coherent sum).
% Stored (both observers): coherent sum E_tot = sum_j E_j and incoherent sum sum_j |E_j|^2, also for the
% two halves of the bunch (odd / even IDs, 1000 electrons each) for a statistical error estimate.
% Also checked: the largest |x|, |y| of any electron inside the field map (must stay < 2 mm).
% Resumable: batches of 200 particles, state saved after each batch.
%
% Syntax: run_case(case_dir, gdf_name, out_file)
%   case_dir : GPT case folder (with initial_distribution.txt and gdf_name), out_file: result .mat

this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
% (read_gpt_particles and particle_traj are in magneto/)                                       % read_gpt_particles, particle_traj
gdf_file = fullfile(case_dir, gdf_name);
N = 2000; batch = 200; nb = N/batch;
obs = {[0 0 11], 11*[sin(0.5e-3) 0 cos(0.5e-3)]}; obs_name = {'onaxis', 'th0p5'};
lambda_nm = logspace(log10(15), log10(150), 5000); lambda = lambda_nm*1e-9; Nw = numel(lambda);
c0 = 299792458;
D0 = readmatrix(fullfile(case_dir, 'initial_distribution.txt'), 'NumHeaderLines', 1);

state_file = [out_file(1:end-4) '_state.mat'];
if exist(state_file, 'file')
    S = load(state_file); E = S.E; I = S.I; done_b = S.done_b; xmax = S.xmax;
    fprintf('Continuing %s: %d of %d batches done\n', case_dir, nnz(done_b), nb);
else
    E = zeros(2, 3, 3, Nw); I = zeros(2, 3, 3, Nw);          % (observer, half: all/odd/even, component, w)
    done_b = false(1, nb); xmax = 0;
end

t0 = tic;
for b = find(~done_b)
    ids = (b-1)*batch + 1 : b*batch; tb = tic;
    tr = read_gpt_particles(gdf_file, ids, true);              % no electric field: skip fE
    % first output is at z ~ 0.65 m: compare with the straight-line drift from the start
    zp = tr.z(:,1); xp = D0(ids,1) + D0(ids,4)./D0(ids,6).*(zp - D0(ids,3)); yp = D0(ids,2) + D0(ids,5)./D0(ids,6).*(zp - D0(ids,3));
    d = max(abs([tr.x(:,1) - xp; tr.y(:,1) - yp]));
    if d > 1e-9; error('ID / initial distribution mismatch (%.2e m)', d); end
    tread = toc(tb);
    for k = 1:numel(ids)
        traj = particle_traj(tr, k);
        inmap = traj.r(3,:) > 0.67996 & traj.r(3,:) < 1.32004;
        xmax = max(xmax, max(abs(traj.r(1:2, inmap)), [], 'all'));
        traj = field_region(traj);
        h = 2 + mod(ids(k), 2);                                 % 2: odd IDs, 3: even IDs (index into 'half')
        for q = 1:2
            s = emission_time_spectrum(traj, obs{q}, lambda, 'fB', norm(obs{q})/c0);
            for hh = [1 h]
                E(q, hh, :, :) = squeeze(E(q, hh, :, :)) + s.Ew;
                I(q, hh, :, :) = squeeze(I(q, hh, :, :)) + abs(s.Ew).^2;
            end
        end
    end
    done_b(b) = true;
    save(state_file, 'E', 'I', 'done_b', 'xmax');
    fprintf('  batch %2d / %d: read %.0f s, total %.0f s; %.1f min elapsed; max |x|,|y| in the map so far %.3f mm\n', b, nb, tread, toc(tb), toc(t0)/60, 1e3*xmax);
    clear tr
end
if xmax > 1.95e-3; warning('an electron reached %.3f mm, close to / beyond the 2 mm map', 1e3*xmax); end
res.E = E; res.I = I; res.lambda_nm = lambda_nm; res.omega = 2*pi*c0./lambda; res.obs = obs; res.obs_name = obs_name;
res.halves = {'all', 'odd IDs', 'even IDs'}; res.xmax_in_map = xmax; res.gdf = gdf_file;
res.sample_sig_d = std(D0(:,7))/mean(D0(:,7)); res.sample_std_xp = std(D0(:,4)./D0(:,6)); res.sample_std_yp = std(D0(:,5)./D0(:,6));
res.note = 'E(observer, half, component, omega): coherent sum of E~_j; I: sum of |E~_j|^2; spectra = eps0 c |obs|^2/pi * |E|^2 or * I';
save(out_file, 'res');
delete(state_file);
fprintf('  done in %.1f min\n', toc(t0)/60);
end

function traj = field_region(traj)
k = find(any(traj.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(traj.t), k(end)+1);
for f = {'t', 'G'};             traj.(f{1}) = traj.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; traj.(f{1}) = traj.(f{1})(:, k); end
end
