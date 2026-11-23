%% compare_translation
% Strict translation model for the multi-particle benchmark (on axis, z = 11 m).
% Without space charge, energy spread or divergence, an electron starting at
% (x0, y0, z0) at t = 0 would - if the field were the same at its offset - move as the
% reference electron (started at the origin), shifted transversely and in time:
%   r_j(t) = r_ref(t + z0/(beta c)) + (x0, y0, 0),   beta_j(t) = beta_ref(t + z0/(beta c))
% (in the drift before the undulator it reaches every z earlier by z0/(beta c)).
% Implementation: the reference samples are relabelled t_i -> t_i - z0/(beta c), the
% positions shifted by (x0, y0, 0), beta and fB kept; then emission_time_spectrum with
% the same common observer-time reference t_ref as run_multi_onaxis.m. No straight
% sections have to be added: only the field region contributes (radiation term), it
% is fully covered, and the Fourier sum runs over each particle's own samples.
% The exact viewing geometry of each electron (observer angle, path length) is thus
% included in both the model and the GPT-based multi-particle result; their difference
% is only the electrons' real trajectory differences in the field at their offsets.
%
% Outputs: results/compare_translation.(mat|txt), figures/compare_translation.(png|fig)

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; eps0 = 8.85418781e-12;

M = load(fullfile(res_dir, 'multi_onaxis.mat')); M = M.multi;
C = load(fullfile(res_dir, 'compare_form_factor.mat')); C = C.cmp;
lam = M.lambda_nm; w = M.omega; N = M.params.N; obs = M.params.obs_point_m; t_ref = M.params.t_ref_s;
pref = eps0 * c0 * norm(obs)^2 / pi;

% reference electron (field region only)
ref = read_gpt_trajectory(fullfile(root_dir, 'gpt_sims', 'single_electron', 'E240MeV_und15p_paraxial', 'result.gdf'));
k = find(any(ref.fB ~= 0, 1)); k = max(1, k(1)-1) : min(numel(ref.t), k(end)+1);
for f = {'t', 'G'};             ref.(f{1}) = ref.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; ref.(f{1}) = ref.(f{1})(:, k); end
beta0 = sqrt(1 - 1/ref.G(1)^2);                          % drift velocity before the undulator

D0 = readmatrix(fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E240MeV_und15p_N2000', 'initial_distribution.txt'), 'NumHeaderLines', 1);
mE = matfile(fullfile(res_dir, 'multi_onaxis_Ej.mat'));

Esum = zeros(3, numel(w)); incoh = zeros(3, numel(w));
dev2 = zeros(N, 1);                                      % per particle: sum_w |E_j,GPT - E_j,trans|^2 / sum_w |E_j,trans|^2
tic;
for b0 = 1:200:N
    ids = b0 : min(b0 + 199, N);
    Eg = double(mE.Ej(ids, 1:3, 1:numel(w)));            % GPT-based fields, 200 x 3 x Nw
    for m = 1:numel(ids)
        j = ids(m);
        tr = ref;
        tr.t = ref.t - D0(j, 3) / (beta0 * c0);
        tr.r = ref.r + [D0(j, 1); D0(j, 2); 0];
        s = emission_time_spectrum(tr, obs, lam*1e-9, 'fB', t_ref);
        Esum = Esum + s.Ew; incoh = incoh + abs(s.Ew).^2;
        d = squeeze(Eg(m, :, :)) - s.Ew;
        dev2(j) = sum(abs(d(:)).^2) / sum(abs(s.Ew(:)).^2);
    end
    fprintf('particles %d-%d done, %.1f min\n', ids(1), ids(end), toc/60);
end

S_tr  = pref * sum(abs(Esum).^2, 1);                     % translation model, coherent sum
S_mul = M.d2W_dwdOmega.total;
nw = round(0.01 / mean(diff(log(lam)))); env = @(S) movmax(S, 2*nw + 1);
e_tr = abs(S_mul - S_tr) ./ env(S_tr);
e_dc = abs(S_mul - C.S_discrete) ./ env(C.S_discrete);   % previous comparison (phase-only model)
e_inc = abs(M.incoherent_sum.total - pref*sum(incoh, 1)) ./ env(pref*sum(incoh, 1));
fprintf('per particle rms |E_GPT - E_trans| / |E_trans|: median %.2e, max %.2e\n', median(sqrt(dev2)), max(sqrt(dev2)));
cy = corrcoef(sqrt(dev2), D0(:,2).^2); cx = corrcoef(sqrt(dev2), D0(:,1).^2);
fprintf('  correlation with y0^2 %.2f, with x0^2 %.2f\n', cy(1,2), cx(1,2));
for band = [15 25; 25 40; 40 70; 70 150]'
    sel = lam >= band(1) & lam < band(2);
    fprintf('  %3d-%3d nm: |multi - translation|/env max %.2e (median %.1e) | previous phase-only model max %.2e | incoherent sums max %.1e\n', ...
        band, max(e_tr(sel)), median(e_tr(sel)), max(e_dc(sel)), max(e_inc(sel)));
end

tm.lambda_nm = lam; tm.omega = w; tm.S_translation = S_tr; tm.S_multi = S_mul; tm.incoherent_translation = pref*sum(incoh, 1);
tm.err_env_translation = e_tr; tm.err_env_phase_only = e_dc; tm.dev_rms_per_particle = sqrt(dev2); tm.Esum_translation = Esum;
tm.note = 'translation model: reference trajectory shifted by (x0,y0,0) and in time by -z0/(beta c); d2W/domega dOmega [J/(rad/s)/sr]';
save(fullfile(res_dir, 'compare_translation.mat'), 'tm');
fid = fopen(fullfile(res_dir, 'compare_translation.txt'), 'wt');
fprintf(fid, '%% On axis (0,0,11) m, d2W/domega dOmega [J/(rad/s)/sr]: 1 lambda [nm]  2 omega [rad/s]  3 MAGNETO multi (GPT trajectories)  4 translation model  5 |3-4|/local envelope\n');
fprintf(fid, '%-14.8f %-14.6e %-14.6e %-14.6e %-14.6e\n', [lam; w; S_mul; S_tr; e_tr]);
fclose(fid);

xt = [15 20 30 40 50 70 100 150];
f = figure('Position', [100 100 950 600]);
loglog(lam, e_dc, 'Color', [0.7 0.7 0.7]); hold on; loglog(lam, e_tr, 'r-', 'LineWidth', 1);
grid on; xlim([15 150]); set(gca, 'XTick', xt);
legend('single electron \times discrete form factor (phase only)', 'strict translation model (trajectory shifted in space and time)', 'Location', 'best');
xlabel('\lambda [nm]'); ylabel('|multi - model| / local envelope (\pm1% in \lambda)');
title('MAGNETO multi-particle (GPT trajectories) vs translation models, on axis at 11 m');
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'compare_translation.png'), 'Resolution', 200); savefig(f, fullfile(fig_dir, 'compare_translation.fig'));
