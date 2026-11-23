function spec = emission_time_spectrum(traj, obs_point, lambda, bdot_method, t_ref)
%% emission_time_spectrum
% Radiation spectrum at one fixed observer point from a sampled trajectory.
% Same physics as MAGNETO (radiation term of the Lienard-Wiechert field with
% the exact n(t), R(t) of every sample), but the Fourier transform is done
% directly in emission time instead of resampling + FFT in observer time:
%
%   E~(w) = sum_i E_i (1 - n_i.b_i) w_i exp(-i w (tobs_i - tobs_1))
%
% which is the exact change of variables dt_obs = (1 - n.b) dt for a fixed
% observer. w_i are trapezoid weights in t. No sampling period, no aliasing,
% any wavelength grid.
%
% Syntax:
%   spec = emission_time_spectrum(traj, obs_point, lambda, bdot_method)
%
% Inputs:
%   traj         - Struct from read_gpt_trajectory (t, r, b, G, fE, fB, q, m, nmacro).
%   obs_point    - Observer position [x, y, z] [m].
%   lambda       - Wavelengths [m] (vector).
%   bdot_method  - How d(beta)/dt is obtained:
%                  'fB'      - Lorentz force from the GPT fields at the same sample (default)
%                  'central' - Central difference (b(i+1)-b(i-1))/(t(i+1)-t(i-1))
%                  'forward' - Forward difference (b(i+1)-b(i))/dt with r, b at i (as in MAGNETO)
%   t_ref        - (optional) common observer-time reference [s]: the phase is
%                  exp(-i w (tobs_i - t_ref)). Needed when the fields of several
%                  particles are summed coherently. Default: tobs_1 of this particle
%                  (only |E~| is then meaningful).
%
% Output (struct):
%   spec.lambda, spec.omega      - Wavelength [m] and angular frequency [rad/s], 1 x Nw
%   spec.Ew                      - Fourier-transformed field [V s/m], 3 x Nw (x; y; z)
%   spec.d2W_dwdOmega            - Energy [J/(rad/s)/sr], struct .total/.x/.y/.z
%   spec.R_obs                   - |obs_point| [m], used to normalise per steradian

if nargin < 4 || isempty(bdot_method); bdot_method = 'fB'; end

% --- Constants (same as MAGNETO) ---
e0 = 1.60217663e-19; c0 = 299792458; eps0 = 8.85418781e-12;
q  = traj.q * traj.nmacro;   % charge of the (macro)particle [C], negative for electrons

t = traj.t; r = traj.r; b = traj.b; G = traj.G;
N = numel(t);
w_trap = [diff(t) 0]/2 + [0 diff(t)]/2;   % trapezoid weights in t

%% 1. d(beta)/dt
switch lower(bdot_method)
    case 'fb'
        % m d(gamma b c)/dt = q (E + c b x B)  ->  db/dt = q/(gamma m c) (E + c b x B - b (b.E))
        k = 1:N;
        bdot = (traj.q ./ (G * traj.m * c0)) .* ...
            (traj.fE + c0 * cross(b, traj.fB, 1) - b .* sum(b .* traj.fE, 1));
        wt = w_trap(k);
    case 'central'
        k = 2:N-1;
        bdot = (b(:, 3:end) - b(:, 1:end-2)) ./ (t(3:end) - t(1:end-2));
        wt = w_trap(k);
    case 'forward'
        k = 1:N-1;
        bdot = diff(b, 1, 2) ./ diff(t);
        wt = diff(t);
    otherwise
        error('Unknown bdot_method: %s', bdot_method);
end

%% 2. Radiation field at the observer (time domain, at each emission sample)
nvec  = obs_point(:) - r(:, k);
R     = sqrt(sum(nvec.^2, 1));
n_hat = nvec ./ R;
one_minus_nb = 1 - sum(n_hat .* b(:, k), 1);
E = (q / (4*pi*eps0)) * cross(n_hat, cross(n_hat - b(:, k), bdot, 1), 1) ./ (c0 * R .* one_minus_nb.^3);
t_obs = t(k) + R / c0;

%% 3. Fourier transform in emission time
omega = 2*pi*c0 ./ lambda(:)';
% subtract a reference time so that the phase keeps full double precision
if nargin < 5 || isempty(t_ref); t_ref = t_obs(1); end
tau = (t_obs - t_ref)';
Ewt = E .* (one_minus_nb .* wt);
Ew  = zeros(3, numel(omega));
chunk = max(1, floor(1e6 / numel(tau)));                % frequencies per block: <= 1e6 complex (~16 MB, ~40 MB peak)
for j0 = 1:chunk:numel(omega)
    jj = j0:min(j0 + chunk - 1, numel(omega));
    Ew(:, jj) = Ewt * exp(-1i * tau * omega(jj));       % 3 x numel(jj)
end

%% 4. Spectral energy density per solid angle
R_obs = norm(obs_point);
pref = eps0 * c0 * R_obs^2 / pi;                        % d2W/dwdOmega = eps0 c R^2 |E~|^2 / pi
comp = pref * abs(Ew).^2;

spec.lambda = lambda(:)';
spec.omega  = omega;
spec.Ew     = Ew;
spec.d2W_dwdOmega.x     = comp(1, :);
spec.d2W_dwdOmega.y     = comp(2, :);
spec.d2W_dwdOmega.z     = comp(3, :);
spec.d2W_dwdOmega.total = sum(comp, 1);
spec.R_obs       = R_obs;
spec.obs_point   = obs_point(:)';
spec.bdot_method = bdot_method;
end
