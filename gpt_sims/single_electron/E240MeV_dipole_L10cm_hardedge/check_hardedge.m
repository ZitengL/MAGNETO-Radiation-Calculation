%% check_hardedge
% Checks of the hard-edge dipole run (b1 = 1e5 1/m) and comparison with b1 = 1e6 1/m (b1_check/)
% and with the Enge-fringe run (../E240MeV_dipole_L10cm, b1 = 200 1/m):
% field profile, bend angle vs asin(e B0 L / p), bending radius, energy, the arc-midpoint tangent
% (direction = half the bend angle) and its intersection with the plane z = 11 m.
clear;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
me_g = 9.10953e-31; qe_g = 1.6021892e-19; c0 = 299792458;    % GPT constants
B0 = 0.65; L = 0.1;
files = {fullfile(this_dir, 'result.gdf'), fullfile(this_dir, 'b1_check', 'result_b1e6.gdf'), ...
         fullfile(fileparts(this_dir), 'E240MeV_dipole_L10cm', 'result.gdf')};
lab = {'hard edge b1 = 1e5', 'hard edge b1 = 1e6', 'Enge b1 = 200 (old)'};
for f = 1:3
    tr = read_gpt_trajectory(files{f});
    z = tr.r(3,:); y = tr.r(2,:); b = tr.b; B = tr.fB(1,:); G = tr.G(1);
    nz = find(B ~= 0); GB = sqrt(G^2 - 1);
    th = atan2(b(2,end), b(3,end)); th_ex = asin(qe_g*B0*L/(GB*me_g*c0));
    % midpoint tangent: where the direction is half the exit angle
    ang = atan2(b(2,:), b(3,:)); k = find(ang >= th/2, 1);
    s = (th/2 - ang(k-1)) / (ang(k) - ang(k-1)); P = tr.r(:, k-1) + s*(tr.r(:, k) - tr.r(:, k-1));
    yobs = P(2) + (11 - P(3))*tan(th/2);
    fprintf('\n%s: %d points, field nonzero z = %.6f .. %.6f m, field at z = 1 m %.6f T, max |fBx| %.6f T\n', lab{f}, numel(z), z(nz(1)), z(nz(end)), interp1(z, B, 1.0), max(abs(B)));
    fprintf('  gamma %.8f (const %d); exit angle %.6e rad, asin(e B0 L/p) = %.6e rad (rel. diff %.1e); rho = p/(e B0) = %.6f m\n', ...
        G, max(tr.G) - min(tr.G) < 1e-9, th, th_ex, th/th_ex - 1, GB*me_g*c0/(qe_g*B0));
    fprintf('  exit at z = %.4f m: y = %.6f mm; max |x| %.1e m\n', z(end), 1e3*y(end), max(abs(tr.r(1,:))));
    fprintf('  midpoint tangent: point (y, z) = (%.6f mm, %.6f m), direction %.6f mrad -> at z = 11 m: y = %.6f m, distance %.6f m\n', ...
        1e3*P(2), P(3), 1e3*th/2, yobs, norm([0; yobs; 11] - P));
end
