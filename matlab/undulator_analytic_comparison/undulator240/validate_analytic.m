%% validate_analytic
% Checks of the implementation undulator_analytic.m (not of the physics):
%   1. A_k(K) of (26.61) against Wiedemann Table 26.1
%   2. (26.46) against a direct numerical time integration of the same model, eq. (26.34)
%      (same approximations, no Bessel expansion) -> tests the Bessel sums and their truncation
%   3. on axis: (26.46) at omega = k omega_1 against the peak value from (26.59)/(26.61),
%      and (26.51) against (26.46) at the peaks
clear;
c0 = 299792458; e0 = 1.6021892e-19; eps0 = 8.85418781e-12;
gam = 470.66442480; K = 0.980309; Np = 15; lamu = 0.030472767;
thu = atan(11*sin(0.5e-3) / (11*cos(0.5e-3) - 1));

%% 1. Table 26.1
for Kt = [1.0 1.4 2.0]
    t = undulator_analytic(100e-9, gam, Kt, Np, lamu, 0, 0);
    fprintf('K = %.1f: A_1, A_3, A_5, A_7 = %s   (book: K=1.0 0.368 0.179 0.055 0.015; 1.4 0.371 0.354 0.219 0.118; 2.0 0.290 0.423 0.413 0.354)\n', Kt, mat2str(t.A_k, 3));
end

%% 2. direct integration of (26.34)
lam = logspace(log10(15), log10(150), 600) * 1e-9; w = 2*pi*c0 ./ lam;
for th = [0 thu]
    a = undulator_analytic(lam, gam, K, Np, lamu, th, 0);
    bbar = a.params.beta_bar; wp = 2*pi/lamu*bbar*c0; D = 1 + K^2/2 + gam^2*th^2; w1 = wp*2*gam^2/D;
    C = 2*K*bbar*gam*th/D; S = K^2*bbar/(4*D);
    t = linspace(-pi*Np/wp, pi*Np/wp, 4000*Np + 1); psi = wp*t;
    Wd = zeros(size(w));
    for j = 1:numel(w)
        X = -1i*(w(j)/w1)*(psi - C*cos(psi) - S*sin(2*psi));
        Ix = trapz(t, (th + K/gam*sin(psi)) .* exp(X));
        Wd(j) = e0^2*bbar^2*w(j)^2/(16*pi^3*eps0*c0) * abs(Ix)^2;
    end
    nw = 3; env = movmax(Wd, 2*nw + 1);
    fprintf('theta = %.4f mrad: |(26.46) - direct (26.34)| / envelope max %.2e (Bessel orders |n| <= %d, u_max %.2f, v_max %.2f)\n', ...
        th*1e3, max(abs(a.W46.total - Wd) ./ env), a.nmax, a.u_max, a.v_max);
end

%% 3. on-axis peaks
a = undulator_analytic(2*pi*c0 ./ ((1:2:7) * undulator_analytic(100e-9, gam, K, Np, lamu, 0, 0).omega1), gam, K, Np, lamu, 0, 0);
fprintf('on axis at omega = k omega_1, k = 1 3 5 7:\n  (26.46)      %s\n  (26.59/61)   %s\n  (26.51) sum  %s\n', ...
    mat2str(a.W46.total, 5), mat2str(a.peak61, 5), mat2str(a.W51.total, 5));
fprintf('  ratios (26.46)/(26.61) - 1: %s;  (26.51)/(26.46) - 1: %s\n', mat2str(a.W46.total ./ a.peak61 - 1, 3), mat2str(a.W51.total ./ a.W46.total - 1, 3));
