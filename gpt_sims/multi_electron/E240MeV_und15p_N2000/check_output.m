%% check_output
% Checks of the 2000-electron GPT run (result.gdf, stats.gdf):
%   - particle count, beam sizes, energy over the run (stats.gdf)
%   - for a batch of 200 particles (IDs 1-200): energy, and how much each
%     trajectory differs from a pure translation of the bunch-average one
%     (natural focusing in y, field non-uniformity in x, timing)
% Output: printed numbers, check_output.png

clear; close all;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(root_dir);                                                % load_gdf
addpath(fullfile(root_dir, 'matlab', 'radiation_multi_electron'));
c0 = 299792458;

%% 1. stats.gdf
S = load_gdf(fullfile(this_dir, 'stats.gdf')); d = S(end).d;
fprintf('stats: %d output times, t = %.6e .. %.6e s, avgz = %.5f .. %.5f m\n', numel(d.time), d.time(1), d.time(end), d.avgz(1), d.avgz(end));
fprintf('numpar: %d .. %d\n', min(d.numpar), max(d.numpar));
fprintf('stdx: %.4f .. %.4f um; stdy: %.4f .. %.4f um; stdz: %.5f .. %.5f nm\n', 1e6*[min(d.stdx) max(d.stdx) min(d.stdy) max(d.stdy)], 1e9*[min(d.stdz) max(d.stdz)]);
fprintf('stdy start / end: %.4f / %.4f um\n', 1e6*d.stdy(1), 1e6*d.stdy(end));
fprintf('avgG: %.8f .. %.8f; max stdG %.2e\n', min(d.avgG), max(d.avgG), max(d.stdG));

%% 2. batch of 200 particles
tic; tr = read_gpt_particles(fullfile(this_dir, 'result.gdf'), 1:200); fprintf('batch read in %.1f s, %d groups\n', toc, numel(tr.t));
fprintf('G: %.10f .. %.10f\n', min(tr.G(:)), max(tr.G(:)));
x0 = tr.x(:, 1); y0 = tr.y(:, 1); z0 = tr.z(:, 1);
dx = tr.x - x0; dy = tr.y - y0; dz = tr.z - z0;              % motion relative to the start of the output
xr = mean(dx, 1); zr = mean(dz, 1);                           % bunch-average (reference) motion
ex = dx - xr; ey = dy; ez = dz - zr;                           % deviation from a pure translation
fprintf('horizontal wiggle amplitude (reference): %.3f um\n', 1e6*max(abs(xr)));
fprintf('max |x deviation| %.3e m; max |y motion| %.3e m (max |y0| %.1f um); max |z deviation| %.3e m\n', ...
    max(abs(ex(:))), max(abs(ey(:))), 1e6*max(abs(y0)), max(abs(ez(:))));
p = polyfit(y0, tr.By(:, end), 1);
fprintf('exit vertical angle vs y0: slope %.4e rad/m (fit), i.e. %.3e rad at y0 = 40 um\n', p(1), p(1)*40e-6);
K = 0.98; ku = 2*pi/0.0304728; G0 = tr.G(1,1); kb = K*ku/(sqrt(2)*G0); Lu = 15*0.0304728;
fprintf('analytic natural focusing: k_beta = %.4f 1/m, -k_beta^2 L_u = %.4e rad/m (15 periods)\n', kb, -kb^2*Lu);
% timing: arrival-time spread through the undulator relative to the initial z offsets
fprintf('max |z deviation| / sigma_z = %.2e\n', max(abs(ez(:))) / 10e-9);

f = figure('Visible', 'off', 'Position', [100 100 900 700]);
z = tr.z(1, :);
subplot(3,1,1); plot(z, 1e6*xr); ylabel('x_{ref} - x_0 [\mum]'); grid on; title('bunch-average horizontal motion');
subplot(3,1,2); plot(z, 1e9*ey(1:20, :)); ylabel('y - y_0 [nm]'); grid on; title('vertical motion of 20 particles (natural focusing)');
subplot(3,1,3); plot(z, 1e12*ex(1:20, :)); ylabel('x deviation [pm]'); xlabel('z [m]'); grid on; title('deviation from a pure translation in x');
exportgraphics(f, fullfile(this_dir, 'check_output.png'), 'Resolution', 150);
