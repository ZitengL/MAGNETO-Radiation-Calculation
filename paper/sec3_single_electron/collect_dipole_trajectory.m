%% collect_dipole_trajectory
% Single-electron trajectory through the short dipole (240 MeV, 0.1 m, 0.65 T, Enge fringe,
% centre z = 1 m), for the geometry panel of the dipole figure. Read from the GPT run used for
% the on-axis spectrum; nothing is recomputed.
% Source (read only): gpt_sims/single_electron/E240MeV_dipole_L10cm/result.gdf (tout 3.75 um)
% Output: data/dipole_trajectory.txt, data/dipole_trajectory.mat
%   columns: z [m], y [m], bend angle atan(by/bz) [rad], |B| at the particle [T]
% Run: matlab -batch "collect_dipole_trajectory"

this_dir = fileparts(mfilename('fullpath'));
root     = fullfile(this_dir, '..', '..');
addpath(fullfile(root, 'magneto'));
tr = read_gpt_trajectory(fullfile(root, 'gpt_sims', 'single_electron', 'E240MeV_dipole_L10cm', 'result.gdf'));
k  = 1:20:numel(tr.t);                                  % every 75 um is enough for a plot
z  = tr.r(3, k)'; y = tr.r(2, k)'; ang = atan2(tr.b(2, k), tr.b(3, k))'; B = vecnorm(tr.fB(:, k))';
fprintf('z = %.4f .. %.4f m, final bend angle %.4f mrad, max |B| = %.5f T, max |x| = %.2e m\n', ...
    z(1), z(end), 1e3*ang(end), max(B), max(abs(tr.r(1, :))));

% tangent points of the two observers (on axis (0,0,11) m; arc midpoint tangent (0,0.406788,11) m)
obs = [0 0 11; 0 0.406788 11];
% and the part of the trajectory that radiates towards each observer: angle between the direction of
% motion and the line of sight < 1/gamma
zr = zeros(2, 2);
for i = 1:2
    d = obs(i, :)' - tr.r; d = d ./ vecnorm(d);
    a = acos(min(1, sum(d .* tr.b, 1) ./ vecnorm(tr.b)));
    [~, j] = min(a);                                    % direction of motion closest to the line of sight
    tp(i, :) = tr.r(:, j)'; %#ok<SAGROW>
    in = find(a < 1/mean(tr.G)); zr(i, :) = tr.r(3, in([1 end]));
    fprintf('observer %d: closest direction at z = %.4f m (y = %.4f mm, |B| = %.4f T); angle < 1/gamma for z = %.4f .. %.4f m (|B| at the end %.4f T)\n', ...
        i, tp(i, 3), 1e3*tp(i, 2), norm(tr.fB(:, j)), zr(i, 1), zr(i, 2), norm(tr.fB(:, in(end))));
end

save(fullfile(this_dir, 'data', 'dipole_trajectory.mat'), 'z', 'y', 'ang', 'B', 'obs', 'tp', 'zr');
fid = fopen(fullfile(this_dir, 'data', 'dipole_trajectory.txt'), 'wt');
fprintf(fid, '%% Single electron, 240 MeV, short dipole 0.1 m / 0.65 T (Enge fringe, gap 1 cm), centre z = 1 m, bend upwards (+y).\n');
fprintf(fid, '%% Source: gpt_sims/single_electron/E240MeV_dipole_L10cm/result.gdf, every 20th output point (75 um).\n');
fprintf(fid, '%% Part radiating towards the observer (angle to the line of sight < 1/gamma): on axis (0,0,11) m -> z = %.4f .. %.4f m (start of the output);\n%% arc midpoint (0,0.406788,11) m -> z = %.4f .. %.4f m, tangent point z = %.4f m, y = %.4e m\n', zr(1,1), zr(1,2), zr(2,1), zr(2,2), tp(2,3), tp(2,2));
fprintf(fid, '%% Columns: 1 z [m]  2 y [m]  3 bend angle [rad]  4 |B| [T]\n');
fprintf(fid, '%.10e %.10e %.10e %.10e\n', [z y ang B]');
fclose(fid);
