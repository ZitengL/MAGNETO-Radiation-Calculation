%% collect_sec5_fringe_vs_D
% Interference at one fixed wavelength as a function of the undulator-dipole distance D
% (single electron, on axis, fine scan D = 0.500-0.700 m in 1 mm steps). Nothing is recomputed.
%   MAGNETO: |E_x +- E_y|^2/2 (the +-45 deg projections), in d2W/domega dOmega [J/(rad/s)/sr],
%            per sr of the 11 m sphere centred at the origin.
%   Prediction from the geometry alone: the two sources keep their spectra, only the delay of the
%   dipole light changes, so |E_x +- E_y|^2/2 = (S_und + S_dip)/2 +- sqrt(S_und S_dip) cos(phi_0 + omega dtau(D)),
%   with S_und, S_dip the two sources computed alone at each D, phi_0 = arg(E_x E_y*) at D = 0.5 m and dtau(D) the predicted delay
%   (light leaving the dipole entrance together with the electron, exact geometry to the observer).
% Source (read only): matlab/interference_und_dipole/results/interference_scan.mat
% Output: data/fringe_vs_D_<lam>nm_onaxis.txt
% Run: matlab -batch "collect_sec5_fringe_vs_D"

this_dir = fileparts(mfilename('fullpath'));
root     = fullfile(this_dir, '..', '..', '..');
data_dir = fullfile(this_dir, 'data');
lam0 = 70;                                               % [nm], suggested by the user (strong fringes, V ~ 0.98)

s = load(fullfile(root, 'matlab', 'interference_und_dipole', 'results', 'interference_scan.mat'));
A = s.A.onaxis; L = s.A.lambda_nm(:); w = s.A.omega(:);
[~, j] = min(abs(L - lam0));
sel = find(A.D >= 0.5 - 1e-9 & A.D <= 0.7 + 1e-9); [D, o] = sort(A.D(sel)); sel = sel(o); D = D(:);
i0 = sel(abs(D - 0.5) < 1e-9);
Sp = A.Sp45(sel, j); Sm = A.Sm45(sel, j); Sref = (A.S_und(1, j) + A.S_dip(sel, j))/2;   % S_und does not depend on D (one row)
phi0 = atan2(-A.S3(i0, j), A.S2(i0, j));                 % S2 = 2Re(Ex Ey*), S3 = -2Im(Ex Ey*)
dt   = A.dtau_pred(sel) - A.dtau_pred(i0);  dt = dt(:);  % [s], relative to D = 0.5 m
% intensities of the two sources alone at each D (the dipole moves towards the observer, so S_dip
% grows by ~4.5% over the scan); only the phase is predicted
Su = A.S_und(1, j); Sd = A.S_dip(sel, j); Sd = Sd(:);
pred = @(sg) (Su + Sd)/2 + sg*sqrt(Su*Sd).*cos(phi0 + w(j)*dt);
Pp = pred(1); Pm = pred(-1);
amp = sqrt(Su*mean(Sd));
fprintf('lambda = %.4f nm: visibility at D = 0.5 m %.4f; fringe period in D %.2f mm (2 gamma^2 lambda = %.2f mm)\n', ...
    L(j), A.visibility45(i0, j), 1e3*L(j)*1e-9/(299792458*c_dtau_slope(dt, D)), 1e3*2*470.66442480^2*L(j)*1e-9);
fprintf('max |MAGNETO - prediction| / fringe amplitude sqrt(SxSy): +45 %.2e, -45 %.2e\n', max(abs(Sp - Pp))/amp, max(abs(Sm - Pm))/amp);
fprintf('S_dip over the scan: %.4e .. %.4e (relative range %.2e); S_und = %.4e\n', min(Sd), max(Sd), (max(Sd) - min(Sd))/mean(Sd), Su);

fid = fopen(fullfile(data_dir, sprintf('fringe_vs_D_%gnm_onaxis.txt', lam0)), 'w');
fprintf(fid, '%% Single electron, on axis (0,0,11) m, lambda = %.6f nm. d2W/domega dOmega [J/(rad/s)/sr], per sr of the 11 m sphere.\n', L(j));
fprintf(fid, '%% |E_x +- E_y|^2/2 = +-45 deg projections (MAGNETO, undulator + dipole computed together); reference (S_und + S_dip)/2.\n');
fprintf(fid, '%% Prediction: (S_und + S_dip)/2 +- sqrt(S_und S_dip) cos(phi_0 + omega dtau_pred(D)); S_und, S_dip computed alone at each D, phi_0 at D = 0.5 m.\n');
fprintf(fid, '%% Columns: 1 D [m]  2 |Ex+Ey|^2/2  3 |Ex-Ey|^2/2  4 (S_und+S_dip)/2  5 prediction +45  6 prediction -45  7 c*dtau_pred rel. D = 0.5 m [um]\n');
fprintf(fid, '%.6f %.10e %.10e %.10e %.10e %.10e %.10e\n', [D Sp Sm Sref Pp Pm 299792458e6*dt]');
fclose(fid);

function k = c_dtau_slope(dt, D)
% d(tau)/dD [s/m] from a straight-line fit
p = polyfit(D, dt, 1); k = p(1);
end
