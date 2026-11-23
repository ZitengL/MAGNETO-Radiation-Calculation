%% check_direct_sigz
% Check of the DIRECT 5e4-particle scan (run_direct_sigz.m / analyze_direct_sigz.m) against the expected value of the
% SAME 5e4-macro-particle beam from run_large_sigz.m (expected value over the longitudinal random numbers):
%   <pref |sum E|^2> = pref [ |F|^2 (|sum E0|^2 - sum |E0|^2) + sum |E0|^2 ],  |F|^2 = exp(-(w sz/c)^2)
% i.e. expected coherent part + macro-particle noise floor (NOT divided by nmacro: the direct sum has 5e4 particles).
% z = (W_direct - W_expected) / jackknife sd of W_direct. Expect |z| of order 1 (one realisation).
% Output: results/direct_sigz/check_vs_expected.txt
clear;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results');
c0 = 299792458; eps0 = 8.85418781e-12;
R = load(fullfile(res_dir, 'large_sigz.mat')); R = R.R;
A = load(fullfile(res_dir, 'direct_sigz', 'summary.mat'));
lam = R.lambda; om = R.omega; pref = eps0*c0*norm(R.obs)^2/pi;
if max(abs(A.lam - lam)) > 0; error('grid differs'); end
E0 = sum(squeeze(R.E0(:,1,1,:) + R.E0(:,2,1,:)), 1); I = sum(squeeze(R.I(:,3,1,:)), 1);
bands = A.bands; bn = {'D1', 'D2'};
Wb = @(S, b) -trapz(om(lam >= bands(b,1) & lam <= bands(b,2)), S(lam >= bands(b,1) & lam <= bands(b,2)));
fid = fopen(fullfile(res_dir, 'direct_sigz', 'check_vs_expected.txt'), 'w');
fprintf(fid, '%% direct W vs expected value of the same 5e4-particle beam (coherent expected + macro noise floor), J/sr\n');
fprintf(fid, '%% sigma_z_nm | D1: W_direct  W_expected  coherent_exp/noise  z | D2: same\n');
for s = 1:numel(A.sz)
    F2 = exp(-(om*A.sz(s)/c0).^2);
    Sexp = pref*(F2.*(abs(E0).^2 - I) + I); Scoh = pref*F2.*(abs(E0).^2 - I);
    fprintf(fid, '%5.0f', A.sz(s)*1e9);
    for b = 1:2
        We = Wb(Sexp, b); Wc = Wb(Scoh, b);
        fprintf(fid, ' | %s %.4e %.4e %8.3f %6.2f', bn{b}, A.W(b,s), We, Wc/A.Wn(b,s), (A.W(b,s) - We)/A.Wsd(b,s));
    end
    fprintf(fid, '\n');
end
fclose(fid); type(fullfile(res_dir, 'direct_sigz', 'check_vs_expected.txt'));
