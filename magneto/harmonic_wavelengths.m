function h = harmonic_wavelengths(traj, harmonics)
%% harmonic_wavelengths
% On-axis resonant wavelengths of the undulator, estimated from the GPT
% trajectory in two ways:
%   'K'    : lambda_n = lambda_u/(2 n gamma^2) (1 + K^2/2), with
%            K = e B0 lambda_u / (2 pi m c)  (= 0.934 B0[T] lambda_u[cm]),
%            B0 = peak on-axis field seen by the electron (max |fB|).
%   'beta' : lambda_1 = lambda_u (1/<beta_z> - 1), the slippage of the light
%            ahead of the electron per period; <beta_z> = Delta z / (c Delta t)
%            over the same central window used for lambda_u. lambda_n = lambda_1/n.
% The central window runs between the 2nd and the 2nd-to-last zero crossings
% of the main field component, so the distorted end half-periods are excluded.
%
% Syntax:
%   h = harmonic_wavelengths(traj)
%   h = harmonic_wavelengths(traj, harmonics)     % default harmonics = [1 3 5]
%
% Output (struct):
%   h.harmonics, h.lambda_u [m], h.B0 [T], h.K_B, h.beta_z_mean, h.z_window [m]
%   h.lambda_K [m], h.lambda_beta [m]   - 1 x numel(harmonics)

if nargin < 2; harmonics = [1 3 5]; end

e0 = 1.60217663e-19; c0 = 299792458; me = 9.1093837e-31;

z = traj.r(3, :);
G = mean(traj.G);
% main field component = the transverse one with the largest peak
[~, ic] = max(max(abs(traj.fB(1:2, :)), [], 2));
B = traj.fB(ic, :);

% zero crossings of the main component between the first and last pole
% (|B| > Bpk/2), so sign changes of the weak fringe tails are ignored
Bpk = max(abs(B));
poles = find(abs(B) > 0.5*Bpk);
inside = false(1, numel(B) - 1);
inside(poles(1):poles(end)-1) = true;
i0 = find(B(1:end-1) .* B(2:end) < 0 & inside);
zc = z(i0) - B(i0) .* (z(i0+1) - z(i0)) ./ (B(i0+1) - B(i0));
if numel(zc) < 4; error('Fewer than 4 zero crossings found'); end
zw = zc(2:end-1);                       % central window
n_half = numel(zw) - 1;
lambda_u = 2 * (zw(end) - zw(1)) / n_half;

% slippage s(z) = c t - z is smooth and monotonic: interpolate at the window ends
s  = c0 * traj.t - z;
sw = interp1(z, s, zw([1 end]), 'pchip');
dz = zw(end) - zw(1);
beta_z_mean = dz / (dz + diff(sw));
lambda1_beta = lambda_u * (1/beta_z_mean - 1);

K_B = e0 * Bpk * lambda_u / (2*pi*me*c0);
lambda1_K = lambda_u / (2*G^2) * (1 + K_B^2/2);

h.harmonics   = harmonics;
h.lambda_u    = lambda_u;
h.n_half_periods_in_window = n_half;
h.z_window    = zw([1 end]);
h.B0          = Bpk;
h.K_B         = K_B;
h.gamma       = G;
h.beta_z_mean = beta_z_mean;
h.lambda_K    = lambda1_K ./ harmonics;
h.lambda_beta = lambda1_beta ./ harmonics;
h.K_from_beta = sqrt(2 * (lambda1_beta * 2*G^2 / lambda_u - 1));   % K that reproduces lambda1_beta
end
