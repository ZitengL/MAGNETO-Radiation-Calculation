%% check_output
% Checks of the 2000-electron dipole run: particle count, sizes, energy (stats.gdf),
% and for IDs 1-100 the departure of each trajectory from a pure translation of the
% bunch-average one (transverse and longitudinal).
clear;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(root_dir); addpath(fullfile(root_dir, 'matlab', 'radiation_multi_electron'));
S = load_gdf(fullfile(this_dir, 'stats.gdf')); d = S(end).d;
fprintf('stats: %d output times, avgz %.5f .. %.5f m, numpar %d .. %d\n', numel(d.time), d.avgz(1), d.avgz(end), min(d.numpar), max(d.numpar));
fprintf('stdx %.4f .. %.4f um, stdy %.4f .. %.4f um, stdz %.5f .. %.5f nm; avgG %.8f, max stdG %.1e\n', 1e6*[min(d.stdx) max(d.stdx) min(d.stdy) max(d.stdy)], ...
    1e9*[min(d.stdz) max(d.stdz)], mean(d.avgG), max(d.stdG));
fprintf('bunch centre at the end: avgy %.4f mm, avgBy %.5e\n', 1e3*d.avgy(end), d.avgBy(end));
tr = read_gpt_particles(fullfile(this_dir, 'result.gdf'), 1:100, true);
x0 = tr.x(:,1); y0 = tr.y(:,1); z0 = tr.z(:,1);
ex = (tr.x - x0) - mean(tr.x - x0, 1); ey = (tr.y - y0) - mean(tr.y - y0, 1); ez = (tr.z - z0) - mean(tr.z - z0, 1);
fprintf('departure from a pure translation: max |x| %.2e m, max |y| %.2e m, max |z| %.2e m (rms at exit %.2e m)\n', ...
    max(abs(ex(:))), max(abs(ey(:))), max(abs(ez(:))), std(ez(:, end)));
fprintf('G range %.10f .. %.10f; max |x| %.1f um, y range %.3f .. %.3f mm\n', min(tr.G(:)), max(tr.G(:)), 1e6*max(abs(tr.x(:))), 1e3*min(tr.y(:)), 1e3*max(tr.y(:)));
