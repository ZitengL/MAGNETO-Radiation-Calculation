%% check_convergence
% tout convergence of the most demanding case (sd0_en12: en = 12 um, largest angles and offsets):
% GPT output step 60 um (results/sd0_en12.mat) vs 30 um (results/sd0_en12_dz30um.mat), same particles.
% Metric (as in the other convergence tests): |S_60 - S_30| / envelope, envelope = moving max of S_30
% over +-1% in lambda; reported for the coherent and incoherent spectra, both observers, overall and
% within +-10% of each harmonic. Output: results/convergence_sd0_en12.txt

clear;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results');
c0 = 299792458; eps0 = 8.85418781e-12;
A = load(fullfile(res_dir, 'sd0_en12.mat')); A = A.res; B = load(fullfile(res_dir, 'sd0_en12_dz30um.mat')); B = B.res;
lam = A.lambda_nm; obs = A.obs; oname = A.obs_name;
lam1 = [101.8283, 106.44]; harm = {[1 3 5], [1 2 4 5]};
fid = fopen(fullfile(res_dir, 'convergence_sd0_en12.txt'), 'wt');
fprintf(fid, '%% sd0_en12, tout 60 um vs 30 um: max |S_60 - S_30| / envelope(S_30, +-1%% in lambda)\n');
for q = 1:2
    pref = eps0*c0*norm(obs{q})^2/pi;
    S = {pref*squeeze(sum(abs(A.E(q,1,:,:)).^2, 3))', pref*squeeze(sum(abs(B.E(q,1,:,:)).^2, 3))'; ...
         pref*squeeze(sum(A.I(q,1,:,:), 3))',         pref*squeeze(sum(B.I(q,1,:,:), 3))'};
    kind = {'coherent', 'incoherent'};
    for k = 1:2
        env = movmax_lam(lam, S{k,2}, 0.01); d = abs(S{k,1} - S{k,2})./env;
        line = sprintf('%-7s %-10s overall %.2e', oname{q}, kind{k}, max(d));
        for n = harm{q}; w = abs(lam/(lam1(q)/n) - 1) < 0.1; line = [line sprintf('   n=%d: %.2e', n, max(d(w)))]; end %#ok<AGROW>
        fprintf('%s\n', line); fprintf(fid, '%s\n', line);
    end
end
fclose(fid);

function env = movmax_lam(lam, S, f)
env = zeros(size(S));
for i = 1:numel(lam); env(i) = max(S(abs(lam/lam(i) - 1) <= f)); end
end
