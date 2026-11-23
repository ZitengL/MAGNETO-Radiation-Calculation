%% check_dipole
% Checks of the 240 MeV short-dipole run: field profile vs Enge model, bend angle,
% bending radius and critical wavelength, energy, Lorentz-force integration.
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
tr = read_gpt_trajectory(fullfile(this_dir, 'result.gdf'));
z = tr.r(3,:); y = tr.r(2,:); b = tr.b; B = tr.fB(1,:); G = tr.G;
e = 1.60217663e-19; me = 9.1093837e-31; c = 299792458;
B0 = 0.65; b1 = 200; zin = 0.95; zout = 1.05;
Bm = -B0 ./ (1 + exp(b1*(zin - z))) ./ (1 + exp(b1*(z - zout)));
nz = find(B ~= 0);
fprintf('output: %d points, z = %.5f .. %.5f m; field nonzero z = %.4f .. %.4f m\n', numel(z), z(1), z(end), z(nz(1)), z(nz(end)));
fprintf('field at centre %.6f T; max |fBx - Enge| (field region) %.2e T; field integral %.6f T m (hard edge 0.065)\n', ...
    interp1(z, B, 1.0), max(abs(B(nz) - Bm(nz))), trapz(z, B));
GB = sqrt(G(1)^2 - 1); rho = GB*me*c/(e*B0);
th = atan2(b(2,end), b(3,end));
fprintf('gamma %.6f (const: %d); rho = %.4f m; exit angle %.5e rad, asin(e*int(B)/p) = %.5e rad\n', G(1), max(G) - min(G) < 1e-9, rho, th, asin(e*abs(trapz(z, B))/(GB*me*c)));
fprintf('lambda_c = 4 pi rho / (3 gamma^3) = %.2f nm; y at exit (z = %.3f m) %.3f mm; max |x| %.1e\n', 4*pi*rho/(3*G(1)^3)*1e9, z(end), 1e3*y(end), max(abs(tr.r(1,:))));
bd = (-e/(G(1)*me)) * cross(b, tr.fB); bi = b(:,1) + cumtrapz(tr.t, bd, 2);
fprintf('Lorentz-force integral check: max |by - by_int| / by_end = %.2e\n', max(abs(b(2,:) - bi(2,:)))/abs(b(2,end)));
fprintf('t_obs monotonic for (0,0,11): %d\n', all(diff(tr.t + sqrt(tr.r(1,:).^2 + y.^2 + (11 - z).^2)/c) > 0));
