function a = undulator_analytic(lambda, gam, K, Np, lamu, theta, phi, kmax)
%% undulator_analytic
% Single-electron spectrum d2W/(domega dOmega) [J/(rad/s)/sr] of an ideal planar undulator,
% Wiedemann, Particle Accelerator Physics, 4th ed., Sec. 26.2 (SI units):
%   (26.46) full double Bessel sum over the Np periods  -> a.W46 (main result)
%   (26.51) single-harmonic form, summed over k = 1..kmax -> a.W51 (check near the peaks only)
%   (26.59)-(26.61) on-axis peak values of the odd harmonics -> a.peak61 (at omega = k omega_1)
% Ideal trajectory: beta_x = -(K/gamma) sin(omega_p t), exactly Np periods, far field; theta is the
% observation angle from the undulator axis, phi the azimuth (phi = 0: wiggle plane).
% As in the book, the third integral of (26.35), eq. (26.36) (straight drift inside the magnet length),
% is NOT included (user's decision); it vanishes on axis.
% Notes on the printed formulas (checked against each other):
%   - the A1 term of (26.46) is used in the form that follows from (26.43):
%       A1 * sum J_m(u) J_n(v) i^n * i (pi Np/omega_p) [sinc(R+1) - sinc(R-1)]
%   - (26.51) is printed with an extra 1/gamma^2 in the prefactor; with A0 = theta cos(phi),
%     A1 = K/gamma, B0 = theta sin(phi) (p. 907) this makes it gamma^2 too small compared with
%     (26.46) and (26.59)/(26.61). It is used here without that factor.
%   sinc(x) = sin(pi Np x)/(pi Np x).
%
% Syntax:
%   a = undulator_analytic(lambda, gam, K, Np, lamu, theta, phi)
%   a = undulator_analytic(lambda, gam, K, Np, lamu, theta, phi, kmax)   % default kmax = 8
% Inputs: lambda [m] (vector), gam, K (peak), Np, lamu [m], theta [rad], phi [rad]
% Output (struct): W46.x/.y/.total, W51.total and W51.k (kmax x Nw), omega1, lambda1, peak61 (odd k),
%                  k_odd, nmax (Bessel truncation), u_max, v_max

if nargin < 8; kmax = 8; end
e0 = 1.6021892e-19;                     % same charge as the GPT particle (MAGNETO uses traj.q)
c0 = 299792458; eps0 = 8.85418781e-12;
w = 2*pi*c0 ./ lambda(:)';

bbar = sqrt(1 - 1/gam^2) * (1 - K^2/(4*gam^2));          % average beta, (26.25)
wp = 2*pi/lamu * bbar * c0;                              % (26.24)
D = 1 + K^2/2 + gam^2*theta^2;
w1 = wp * 2*gam^2 / D;                                   % (26.31)
C = 2*K*bbar*gam*theta*cos(phi) / D;                     % (26.37a)
S = K^2*bbar / (4*D);                                    % (26.37b)
A0 = theta*cos(phi); A1 = K/gam; B0 = theta*sin(phi);
pref = e0^2 * bbar^2 * w.^2 / (16*pi^3*eps0*c0);         % a = r_c m c beta^2 omega^2/(4 pi^2), SI
x = w / w1; u = x*S; v = x*C;
sincN = @(y) sinc_np(y, Np);

%% (26.46)
nmax = ceil(max([abs(u) abs(v)])) + 25;                  % J_n(z) negligible for |n| >> z
Ax = zeros(size(w)); Ay = zeros(size(w));
L = 2*pi*Np/wp;
idx = -2*nmax-2 : 2*nmax+2;                              % all orders needed below
[NU, ZU] = ndgrid(idx, u); Ju = besselj(NU, ZU);          % (orders) x Nw
[NV, ZV] = ndgrid(idx, v); Jv = besselj(NV, ZV);
J = @(T, n) T(n - idx(1) + 1, :);
for m = -nmax:nmax
    Jm = J(Ju, m);
    for n = -nmax:nmax
        c = Jm .* J(Jv, n) * 1i^n;
        R = x - n - 2*m;
        s0 = L * sincN(R);
        s1 = 1i*(pi*Np/wp) * (sincN(R + 1) - sincN(R - 1));
        Ax = Ax + c .* (A0*s0 + A1*s1);
        Ay = Ay + c .* (B0*s0);
    end
end
a.W46.x = pref .* abs(Ax).^2; a.W46.y = pref .* abs(Ay).^2; a.W46.total = a.W46.x + a.W46.y;

%% (26.51), harmonic by harmonic
a.W51.k = zeros(kmax, numel(w));
mm = -nmax:nmax;
bj = @(n, z) besselj(n, z);                              % orders up to |k| + 2 nmax + 1
for k = 1:kmax
    Sx0 = 0; Sy0 = 0; Sp = 0; Sm = 0;
    for m = mm
        Jm = J(Ju, m);
        Sx0 = Sx0 + Jm .* bj(k - 2*m, v)     * 1i^(k - 2*m);
        Sp  = Sp  + Jm .* bj(k - 2*m + 1, v) * 1i^(k - 2*m + 1);
        Sm  = Sm  + Jm .* bj(k - 2*m - 1, v) * 1i^(k - 2*m - 1);
    end
    Sy0 = Sx0;
    amp_x = A0*Sx0 + 1i*0.5*A1*Sp - 1i*0.5*A1*Sm;
    amp_y = B0*Sy0;
    a.W51.k(k, :) = pref .* L^2 .* sincN(x - k).^2 .* (abs(amp_x).^2 + abs(amp_y).^2);
end
a.W51.total = sum(a.W51.k, 1);

%% (26.59)-(26.61): on-axis peak d2W/domega dOmega at omega = k omega_1 (odd k)
ko = 1:2:kmax;
zeta = ko*K^2 / (4 + 2*K^2);
JJ = besselj((ko-1)/2, zeta) - besselj((ko+1)/2, zeta);   % (26.60)
Ak = ko.^2 * K^2 / (1 + K^2/2)^2 .* JJ.^2;               % (26.61)
a.k_odd = ko; a.A_k = Ak;
a.peak61 = e0^2 * gam^2 * Np^2 / (4*pi*eps0*c0) * Ak;    % hbar*alpha*gamma^2*Np^2*A_k, from (26.59)

a.omega1 = w1; a.lambda1 = 2*pi*c0/w1; a.nmax = nmax; a.u_max = max(abs(u)); a.v_max = max(abs(v));
a.params = struct('gamma', gam, 'K', K, 'Np', Np, 'lambda_u', lamu, 'theta', theta, 'phi', phi, 'beta_bar', bbar);
end

function s = sinc_np(y, Np)
s = ones(size(y)); k = abs(y) > 1e-12;
s(k) = sin(pi*Np*y(k)) ./ (pi*Np*y(k));
end
