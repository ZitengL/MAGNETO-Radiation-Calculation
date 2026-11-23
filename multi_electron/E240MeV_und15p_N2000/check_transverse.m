%% check_transverse
% What drives the departures from a pure translation (check_output.m)?
% Fits the exit horizontal angle against x0 and the longitudinal deviation at
% the exit against x0^2 and y0^2 (field non-uniformity across the beam), using
% the first 200 particles. Also the field each particle sees vs its offset.
clear;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'matlab', 'radiation_multi_electron'));
tr = read_gpt_particles(fullfile(this_dir, 'result.gdf'), 1:200);
x0 = tr.x(:, 1); y0 = tr.y(:, 1); z0 = tr.z(:, 1);
dz = tr.z - z0; ez = dz - mean(dz, 1);
ex = (tr.x - x0) - mean(tr.x - x0, 1);
p = polyfit(x0, tr.Bx(:, end) - mean(tr.Bx(:, end)), 1);
fprintf('exit horizontal angle vs x0: slope %.4e rad/m (fit residual rms %.2e rad)\n', p(1), std(tr.Bx(:, end) - mean(tr.Bx(:, end)) - polyval(p, x0)));
A = [x0.^2, y0.^2, ones(size(x0))];
c = A \ ez(:, end);
fprintf('exit z deviation = %.4e * x0^2 + %.4e * y0^2 + %.2e  [m]; residual rms %.2e m\n', c, std(ez(:, end) - A*c));
fprintf('  -> at x0 = 40 um: %.3f nm; at y0 = 40 um: %.3f nm\n', c(1)*1.6e-9*1e9, c(2)*1.6e-9*1e9);
fprintf('rms of exit z deviation over the 200 particles: %.3f nm (sigma_z = 10 nm)\n', 1e9*std(ez(:, end)));
% peak field seen vs offset
Bpk = max(abs(tr.fBy), [], 2);
q = [x0.^2, y0.^2, ones(size(x0))] \ Bpk;
fprintf('peak |By| seen = %.6f T %+.4e T/m^2 * x0^2 %+.4e T/m^2 * y0^2\n', q(3), q(1), q(2));
fprintf('  -> relative change at 40 um: x %.2e, y %.2e\n', q(1)*1.6e-9/q(3), q(2)*1.6e-9/q(3));
% the near-axis map covers |x|, |y| <= 0.25 mm: bound for all 2000 particles
D0 = readmatrix(fullfile(this_dir, 'initial_distribution.txt'), 'NumHeaderLines', 1);
fprintf('all 2000 particles: max |x0| %.1f um, max |y0| %.1f um; batch: max |x - x0| %.2f um, max |y - y0| %.3f um\n', ...
    1e6*max(abs(D0(:,1))), 1e6*max(abs(D0(:,2))), 1e6*max(abs(tr.x - x0), [], 'all'), 1e6*max(abs(tr.y - y0), [], 'all'));
