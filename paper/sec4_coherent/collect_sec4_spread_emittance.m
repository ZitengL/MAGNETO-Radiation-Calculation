%% collect_sec4_spread_emittance
% Collects the energy-spread / emittance study (2000 electrons, 15-period undulator, 240 MeV,
% sigma_z = 10 nm defined at t = 0 when the bunch centre is at the origin, 1 m upstream of the
% undulator centre) into ./data/spread_emittance, and copies the original figures into
% ./figures_original/spread_emittance. The spectra are not recomputed; new here are
%   - the coherence ratio r = int S_coh domega / int S_inc domega over the natural-line window of each
%     harmonic (the incoherent FWHM window of the baseline case sd0_en0, the same for all cases;
%     integrated, so that it is not biased by the speckle of S_coh when r ~ 1);
%   - the form-factor-model prediction of the same ratio, r_FF = int S_inc R_FF domega / int S_inc domega,
%     with R_FF(omega) = |sum_j exp(-i omega dt_j)|^2 / N computed from the initial distribution
%     (discrete, this sample) and R_FF = 1 + (N-1)|F|^2 (analytic Gaussian, sigma_z = 10 nm);
%     dt_j = -z_j/(beta c) + (|O - r_j| - |O - r_0|)/c, r_j = (x_w, y_w, z_u) the waist position at the
%     undulator centre. The initial z, x_w, y_w are identical in all cases, so R_FF is the same for all
%     cases: the form-factor model does not see the energy spread or the emittance.
%   - flags for the figures (user's decisions, 2026-10-27):
%     show_coh = r >= 10 (coherent FWHM shown only where coherence dominates);
%     show_inc = false for en = 10, 12 um except n = 1 (FWHM not defined there).
% Sources (read only): matlab/radiation_multi_electron/und15p_N2000_spread_emittance/results/{analysis.mat, sd0_en0.mat}
%                      gpt_sims/multi_electron/E240MeV_und15p_N2000_spread_emittance/<case>/initial_distribution.txt
% Run: matlab -batch "collect_sec4_spread_emittance"

this_dir = fileparts(mfilename('fullpath'));
root     = fullfile(this_dir, '..', '..');
src_dir  = fullfile(root, 'matlab', 'radiation_multi_electron', 'und15p_N2000_spread_emittance');
gdir     = fullfile(root, 'gpt_sims', 'multi_electron', 'E240MeV_und15p_N2000_spread_emittance');
data_dir = fullfile(this_dir, 'data', 'spread_emittance');               if ~exist(data_dir, 'dir'); mkdir(data_dir); end
fig_dir  = fullfile(this_dir, 'figures_original', 'spread_emittance');   if ~exist(fig_dir, 'dir');  mkdir(fig_dir);  end

c0 = 299792458; gam0 = 470.66442480; beta0 = sqrt(1 - 1/gam0^2); zu = 1.0;
sig_z = 10e-9; sig_x = 20e-6; r_min = 10;
A  = load(fullfile(src_dir, 'results', 'analysis.mat'));  W = A.W; lam = A.lam(:); harm = A.harm;
R0 = load(fullfile(src_dir, 'results', 'sd0_en0.mat'));   obs = R0.res.obs; oname = R0.res.obs_name;
omega = 2*pi*c0 ./ (lam*1e-9);
assert(max(abs(R0.res.lambda_nm(:) - lam) ./ lam) < 1e-14, 'lambda grids differ');

%% form factor from the initial distribution (identical in all cases; checked)
names = {W.(oname{1}).name}; nc = numel(names);
D0 = readmatrix(fullfile(gdir, 'sd0_en0', 'initial_distribution.txt'), 'NumHeaderLines', 1);
waist = @(D) [D(:,1) + D(:,4)./D(:,6)*zu, D(:,2) + D(:,5)./D(:,6)*zu, D(:,3)];
P0 = waist(D0); N = size(P0, 1);
for i = 1:nc
    D = readmatrix(fullfile(gdir, names{i}, 'initial_distribution.txt'), 'NumHeaderLines', 1);
    dP = max(abs(waist(D) - P0), [], 1);
    assert(all(dP < [1e-12 1e-12 1e-18]), '%s: waist x, y or initial z differ from sd0_en0', names{i});
end
fprintf('initial z and waist x, y identical in all %d cases; N = %d, sample sigma_z = %.4f nm\n', nc, N, 1e9*std(P0(:,3)));

RFF = struct();
for q = 1:2
    O = obs{q}(:)'; r0 = [0 0 zu];
    Lj = vecnorm(O - [P0(:,1:2), zu*ones(N,1)], 2, 2); L0 = norm(O - r0);
    dpath = (Lj.^2 - L0^2) ./ (Lj + L0);                               % |O - r_j| - |O - r_0|, without cancellation
    dt = -P0(:,3)/(beta0*c0) + dpath/c0;
    Fd = zeros(size(omega));
    for k = 1:numel(omega); Fd(k) = abs(sum(exp(-1i*omega(k)*dt)))^2 / N; end
    th_u = atan(O(1)/(O(3) - zu));                                     % observer angle seen from the undulator centre
    F2 = exp(-(omega/c0).^2 * (sig_z^2/beta0^2 + (sig_x*sin(th_u))^2));
    RFF.(oname{q}) = struct('discrete', Fd, 'gauss', 1 + (N - 1)*F2, 'F2', F2, 'th_u', th_u);
end

%% coherence ratio per case and harmonic
T = struct();
for q = 1:2
    Wq = W.(oname{q}); hl = harm{q}; Rq = RFF.(oname{q});
    lam1 = A.lam1(q); i0 = find(strcmp(names, 'sd0_en0'));
    for i = 1:nc
        r = Wq(i); r.S_inc = r.S_inc(:); r.S_coh = r.S_coh(:);
        for n = 1:numel(hl)
            [la, lb] = fwhm_window(lam, Wq(i0).S_inc(:), lam1/hl(n));     % natural line: FWHM window of the baseline case
            w = lam >= la & lam <= lb;
            I = @(f) abs(trapz(omega(w), f(w)));
            T.(oname{q}).r_int(i,n)   = I(r.S_coh) / I(r.S_inc);
            T.(oname{q}).r_FF_d(i,n)  = I(r.S_inc .* Rq.discrete) / I(r.S_inc);
            T.(oname{q}).r_FF_g(i,n)  = I(r.S_inc .* Rq.gauss) / I(r.S_inc);
            T.(oname{q}).r_peak(i,n)  = r.pk_coh(n) / r.pk_inc(n);
            T.(oname{q}).win(i,n,:)   = [la lb];
            T.(oname{q}).show_coh(i,n) = T.(oname{q}).r_int(i,n) >= r_min;
            T.(oname{q}).show_inc(i,n) = ~(r.en >= 10e-6 && hl(n) ~= 1);
        end
    end
end

%% textbook linewidth: quadrature sum of FWHM contributions (on axis only; NaN at 0.5 mrad)
% lambda_n = lambda_u (1 + K^2/2 + gamma^2 theta^2) / (2 n gamma^2), all terms as FWHM / lambda:
%   natural (sinc^2 line, N periods)      0.886/(n N)
%   energy spread (d lambda/lambda = -2 delta, Gaussian)   2 * 2.355 sigma_delta
%   emittance (d lambda/lambda = gamma^2 theta^2/(1 + K^2/2); theta^2 = theta_x^2 + theta_y^2 is
%   exponentially distributed with mean 2 sigma_theta^2, FWHM of that = ln2 * mean)
%                                          2 ln2 gamma^2 sigma_theta^2 / (1 + K^2/2),  sigma_theta = en/(gamma sigma_x)
% FWHM_quad = sqrt(natural^2 + spread^2 + emittance^2) (exact only for Gaussian shapes).
N_u = 15; K_u = 0.980309;
for q = 1:2
    hl = harm{q};
    for i = 1:nc
        r = W.(oname{q})(i); sth = r.en/(gam0*sig_x);
        for n = 1:numel(hl)
            T.(oname{q}).fwhm_nat(i,n)  = 0.886/(hl(n)*N_u);
            T.(oname{q}).fwhm_quad(i,n) = NaN;
            if q == 1
                T.(oname{q}).fwhm_quad(i,n) = sqrt((0.886/(hl(n)*N_u))^2 + (2*2.3548*r.sig_d)^2 + (2*log(2)*gam0^2*sth^2/(1 + K_u^2/2))^2);
            end
        end
    end
end
%% text output
hdr = sprintf(['%% 2000 electrons, 15-period undulator (lambda_u = 30.47 mm, K = 0.98), 240 MeV, sigma_z = 10 nm at t = 0 (bunch centre at the origin,\n' ...
    '%% 1 m upstream of the undulator centre), waist sigma_x = sigma_y = 20 um at the undulator centre z = 1 m. Spectra d2W/domega dOmega [J/(rad/s)/sr],\n' ...
    '%% per sr of the 11 m sphere centred at the origin. Observer %%s.\n']);
fid = fopen(fullfile(data_dir, 'widths_coherence.txt'), 'wt');
fprintf(fid, '%% Energy-spread / emittance study. FWHM relative to the peak wavelength of the incoherent sum.\n');
fprintf(fid, '%% inc = sum_j|E_j|^2, coh = |sum_j E_j|^2, est = simple estimate (each electron shifted to its own resonance, angular intensity change NOT included).\n');
fprintf(fid, '%% err = half the difference of the two halves (odd / even IDs).\n');
fprintf(fid, '%% r_int = int S_coh domega / int S_inc domega over the natural-line window [lam_a, lam_b] (incoherent FWHM of the baseline case, same for all cases); r_peak = peak S_coh / peak S_inc (biased by speckle when r ~ 1).\n');
fprintf(fid, '%% r_FF_discrete, r_FF_gauss = the form-factor-model prediction of r_int (same window, weighted with S_inc): discrete = this sample, gauss = sigma_z = 10 nm.\n');
fprintf(fid, '%% show_inc = 0: FWHM not defined (en >= 10 um, n > 1); show_coh = 1 if r_int >= %g.\n', r_min);
fprintf(fid, '%% FWHM_quad = textbook quadrature sum sqrt[(0.886/(nN))^2 + (2*2.355 sigma_delta)^2 + (2 ln2 gamma^2 sigma_theta^2/(1+K^2/2))^2], N = 15, on axis only (NaN at 0.5 mrad).\n');
for q = 1:2
    Wq = W.(oname{q}); hl = harm{q}; Tq = T.(oname{q});
    fprintf(fid, '\n%% observer %s (theta_x = %.2f mrad at the origin, %.4f mrad seen from the undulator centre)\n', oname{q}, 1e3*atan(obs{q}(1)/obs{q}(3)), 1e3*RFF.(oname{q}).th_u);
    fprintf(fid, '%% case         sig_d    en_um  n  lam_peak_nm   FWHM_inc   err_inc    FWHM_coh   err_coh    FWHM_est   lam_a_nm    lam_b_nm    r_int       r_peak      r_FF_discrete r_FF_gauss  show_inc show_coh  FWHM_quad\n');
    for i = 1:nc
        r = Wq(i);
        for n = 1:numel(hl)
            fprintf(fid, '%-12s  %-7.4g  %-5.3g  %d  %11.5f  %9.3e  %9.2e  %9.3e  %9.2e  %9.3e  %10.5f  %10.5f  %10.4e  %10.4e  %10.4e    %10.4e  %d        %d         %9.3e\n', ...
                r.name, r.sig_d, 1e6*r.en, hl(n), r.lpk_inc(n), r.fwhm_inc(n), r.err_inc(n), r.fwhm_coh(n), r.err_coh(n), r.fwhm_est(n), ...
                Tq.win(i,n,1), Tq.win(i,n,2), Tq.r_int(i,n), Tq.r_peak(i,n), Tq.r_FF_d(i,n), Tq.r_FF_g(i,n), Tq.show_inc(i,n), Tq.show_coh(i,n), Tq.fwhm_quad(i,n));
        end
    end
end
fclose(fid);

kinds = {'inc', 'S_inc', 'incoherent sum sum_j|E_j|^2'; 'coh', 'S_coh', 'coherent |sum_j E_j|^2'; 'est', 'S_est', 'simple estimate (baseline incoherent line shape shifted per electron)'};
for q = 1:2
    Wq = W.(oname{q});
    for k = 1:size(kinds, 1)
        M = [lam, omega, cell2mat(arrayfun(@(r) r.(kinds{k,2})(:), Wq, 'UniformOutput', false))];
        fid = fopen(fullfile(data_dir, sprintf('spectra_%s_%s.txt', kinds{k,1}, oname{q})), 'wt');
        fprintf(fid, hdr, oname{q});
        fprintf(fid, '%% %s. Columns: 1 lambda_nm  2 omega_rad_s  then one column per case: %s\n', kinds{k,3}, strjoin(names, '  '));
        fprintf(fid, [repmat('%.10e ', 1, size(M, 2) - 1) '%.10e\n'], M');
        fclose(fid);
    end
    Rq = RFF.(oname{q});
    fid = fopen(fullfile(data_dir, sprintf('form_factor_%s.txt', oname{q})), 'wt');
    fprintf(fid, hdr, oname{q});
    fprintf(fid, '%% Form-factor model, ratio R_FF = S_coh/S_inc. discrete: |sum_j exp(-i omega dt_j)|^2/N of this sample; gauss: 1 + (N-1)|F|^2,\n');
    fprintf(fid, '%% |F|^2 = exp[-(omega/c)^2 (sigma_z^2/beta^2 + sigma_x^2 sin^2 theta_u)], sigma_z = 10 nm, sigma_x = 20 um, theta_u = %.4f mrad.\n', 1e3*Rq.th_u);
    fprintf(fid, '%% Columns: 1 lambda_nm  2 omega_rad_s  3 R_FF_discrete  4 R_FF_gauss  5 |F|^2_gauss\n');
    fprintf(fid, '%.10e %.10e %.10e %.10e %.10e\n', [lam, omega, Rq.discrete, Rq.gauss, Rq.F2]');
    fclose(fid);
end

cases_info = struct('name', names, 'sig_d', {W.(oname{1}).sig_d}, 'en', {W.(oname{1}).en});
save(fullfile(data_dir, 'sec4_spread_emittance.mat'), 'W', 'T', 'RFF', 'lam', 'omega', 'harm', 'obs', 'oname', 'cases_info', 'N', 'r_min', 'sig_z', 'sig_x');

%% copy original figures
pngs = dir(fullfile(src_dir, 'figures', '*.png'));
for k = 1:numel(pngs); copyfile(fullfile(pngs(k).folder, pngs(k).name), fullfile(fig_dir, pngs(k).name)); end
copyfile(fullfile(src_dir, 'results', 'convergence_sd0_en12.txt'), fullfile(data_dir, 'convergence_sd0_en12.txt'));

%% console summary
for q = 1:2
    hl = harm{q}; Tq = T.(oname{q});
    fprintf('\n=== %s: r_int (r_FF_discrete, r_FF_gauss)\n', oname{q});
    for i = 1:nc
        fprintf('%-12s', names{i});
        for n = 1:numel(hl); fprintf('  n=%d: %9.3g (%7.3g, %7.3g)', hl(n), Tq.r_int(i,n), Tq.r_FF_d(i,n), Tq.r_FF_g(i,n)); end
        fprintf('\n');
    end
end

function [la, lb] = fwhm_window(lam, S, lc)
% same peak search and half-maximum crossings as analyze_cases.m (main lobe of the peak within +-10% of lc)
win = find(abs(lam/lc - 1) < 0.1); [~, j] = max(S(win)); j = win(j);
x = log(lam(j-1:j+1)); y = S(j-1:j+1); p = polyfit(x - x(2), y, 2); dx = -p(2)/(2*p(1));
if abs(dx) > abs(x(3) - x(2)); dx = 0; end
h = polyval(p, dx)/2;
a = j; while a > 1 && S(a) > h; a = a - 1; end
b = j; while b < numel(S) && S(b) > h; b = b + 1; end
la = interp1(S([a a+1]), lam([a a+1]), h); lb = interp1(S([b-1 b]), lam([b-1 b]), h);
end
