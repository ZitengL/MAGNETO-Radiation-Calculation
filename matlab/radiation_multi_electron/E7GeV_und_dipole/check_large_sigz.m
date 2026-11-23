%% check_large_sigz
% Extra checks of run_large_sigz.m / analyze_large_sigz.m (2026-11-05):
%  (a) which part makes the incoherent level (undulator / dipole self terms, x and y components)
%  (b) check 1 per frequency: residual r = (|Er|^2 - <|E|^2>) / sd, sd^2 = 2 |A|^2 s + s^2,
%      |A|^2 = |F|^2 (|E0|^2 - I), s = I (1 - |F|^2)  (complex Gaussian random part; CLT with 5e4 macro-particles).
%      If the expected-value formula is right, mean(r) ~ 0 and rms(r) ~ 1 (neighbouring frequencies are correlated).
%  (c) where the lambda-grid check (every 2nd point) is largest
%  (d) jackknife sd of the sigma_z where coherent = 10, 1, 0.1 x incoherent
% Output: results/check_large_sigz.txt
clear;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results');
c0 = 299792458; eps0 = 8.85418781e-12;
R = load(fullfile(res_dir, 'large_sigz.mat')); R = R.R;
A = load(fullfile(res_dir, 'large_sigz_summary.mat'));
lam = R.lambda; om = R.omega; k = om/c0; pref = eps0*c0*norm(R.obs)^2/pi; nm = R.nmacro; nch = numel(R.done);
bands = [50e-9 280e-9; 300e-9 1000e-9]; bn = {'D1 50-280 nm', 'D2 300-1000 nm'};
Wb = @(S, b) -trapz(om(lam >= bands(b,1) & lam <= bands(b,2)), S(lam >= bands(b,1) & lam <= bands(b,2)));
fid = fopen(fullfile(res_dir, 'check_large_sigz.txt'), 'w');

% (a) incoherent level of the real beam, sum_j |E_j|^2 / nmacro, by part and component
fprintf(fid, '%% (a) incoherent level (real beam) [J/sr]: part (und / dip / und+dip) x component (x / y)\n');
pn = {'undulator', 'dipole', 'und + dip'}; cn = 'xyz';
for p = 1:3
    for cc = 1:2
        I = squeeze(sum(R.I(:,p,cc,:), 1))'*pref/nm;
        fprintf(fid, '%%   %-9s %s: D1 %.3e, D2 %.3e\n', pn{p}, cn(cc), Wb(I, 1), Wb(I, 2));
    end
end

% (b) per-frequency residual of the realisation
E0 = squeeze(sum(R.E0(:,1,1,:) + R.E0(:,2,1,:), 1)).'; I = squeeze(sum(R.I(:,3,1,:), 1)).';
fprintf(fid, '%% (b) realisation vs expected value per frequency, r = (|Er|^2 - E) / sd: mean(r), rms(r) in each band\n');
for s = 1:numel(R.sz_real)
    F2 = exp(-(k*R.sz_real(s)).^2); A2 = F2.*(abs(E0).^2 - I); sv = I.*(1 - F2);
    Er = squeeze(sum(R.Er(:,s,1,1,:) + R.Er(:,s,2,1,:), 1)).';
    r = (abs(Er).^2 - (A2 + I))./sqrt(2*max(A2, 0).*sv + sv.^2);
    fprintf(fid, '%%   sz %3.0f nm', R.sz_real(s)*1e9);
    for b = 1:2
        m = lam >= bands(b,1) & lam <= bands(b,2);
        fprintf(fid, ' | %s: mean %+.2f, rms %.2f', bn{b}, mean(r(m)), sqrt(mean(r(m).^2)));
    end
    fprintf(fid, '\n');
end

% (c) lambda grid
for b = 1:2
    [g, i] = max(abs(A.gridchk(b,:)));
    fprintf(fid, '%% (c) %s: largest grid change %.1e at sigma_z = %.0f nm; at 10 nm %.1e, at 1000 nm %.1e\n', ...
        bn{b}, g, A.sz(i)*1e9, abs(A.gridchk(b,1)), abs(A.gridchk(b,end)));
end

% (d) jackknife of the threshold sigma_z (linear interpolation of log(coh/inc) on the 10 nm grid)
E0c = squeeze(R.E0(:,1,1,:) + R.E0(:,2,1,:)); Ic = squeeze(R.I(:,3,1,:)); sz = A.sz; f = nch/(nch-1);
thr = @(E0, I, nmac, b, r) cross(arrayfun(@(s) Wb(pref*exp(-(k*s).^2).*(abs(E0).^2 - I), b), sz)/Wb(pref*I/nmac, b), r, sz);
fprintf(fid, '%% (d) sigma_z [nm] where coherent = r x incoherent (interpolated), full sample and jackknife sd\n');
for b = 1:2
    for r = [10 1 0.1]
        full = thr(sum(E0c, 1), sum(Ic, 1), nm, b, r);
        tj = zeros(1, nch);
        for c = 1:nch
            m = setdiff(1:nch, c); tj(c) = thr(f*sum(E0c(m,:), 1), f^2*sum(Ic(m,:), 1), nm*f, b, r);
        end
        fprintf(fid, '%%   %s, r = %4g: %.0f +- %.0f nm\n', bn{b}, r, full*1e9, 1e9*sqrt((nch-1)/nch*sum((tj - mean(tj)).^2)));
    end
end
fclose(fid); type(fullfile(res_dir, 'check_large_sigz.txt'));

function x = cross(q, r, sz)
% first sigma_z where q falls below r (log interpolation)
i = find(q < r, 1); if isempty(i) || i == 1; x = NaN; return; end
x = interp1(log([q(i-1) q(i)]), sz([i-1 i]), log(r));
end
