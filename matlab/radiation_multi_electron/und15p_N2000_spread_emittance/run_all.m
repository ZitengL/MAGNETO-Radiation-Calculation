%% run_all
% Runs all cases of the energy-spread / emittance study (gpt_sims/multi_electron/
% E240MeV_und15p_N2000_spread_emittance/cases.m): GPT (run.bat of the case) if result.gdf is missing,
% then MAGNETO (run_case.m). Skips cases whose results/<case>.mat exists. Finally the tout convergence
% run of the most demanding case (en = 12 um, largest angles): output step 30 um instead of 60 um.
% Disk: before each GPT run, if less than 10 GB are free, the result.gdf of an already processed
% case is deleted (allowed by the user when the disk is short); this is logged.

clear;
this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(fileparts(this_dir)));
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E240MeV_und15p_N2000_spread_emittance');
addpath(gdir);
res_dir = fullfile(this_dir, 'results'); if ~exist(res_dir, 'dir'); mkdir(res_dir); end
C = cases();
jobs = struct('dir', {}, 'gdf', {}, 'out', {}, 'dist', {});
for c = C; jobs(end+1) = struct('dir', fullfile(gdir, c.name), 'gdf', 'result.gdf', 'out', fullfile(res_dir, [c.name '.mat']), 'dist', ''); end %#ok<SAGROW>
jobs(end+1) = struct('dir', fullfile(gdir, 'sd0_en12', 'tout_convergence'), 'gdf', 'result_dz30um.gdf', ...
                     'out', fullfile(res_dir, 'sd0_en12_dz30um.mat'), 'dist', fullfile(gdir, 'sd0_en12'));
t0 = tic;
for i = 1:numel(jobs)
    J = jobs(i);
    if exist(J.out, 'file'); continue; end
    fprintf('\n=== %s (%.1f h elapsed)\n', J.out, toc(t0)/3600);
    if ~exist(fullfile(J.dir, J.gdf), 'file')
        free_disk(gdir, res_dir, C);
        bat = fullfile(J.dir, 'run.bat'); tg = tic;
        [st, msg] = system(sprintf('cmd /c "%s"', bat)); if st ~= 0; error('GPT failed in %s: %s', J.dir, msg); end
        fprintf('  GPT %.1f min\n', toc(tg)/60);
    end
    if ~isempty(J.dist)                           % convergence copy uses the distribution of its parent case
        if ~exist(fullfile(J.dir, 'initial_distribution.txt'), 'file'); copyfile(fullfile(J.dist, 'initial_distribution.txt'), J.dir); end
    end
    run_case(J.dir, J.gdf, J.out);
end
fprintf('\nAll done in %.1f h\n', toc(t0)/3600);

function free_disk(gdir, res_dir, C)
f = java.io.File(gdir); gb = f.getFreeSpace/1e9;
fprintf('  free disk %.1f GB\n', gb);
if gb >= 10; return; end
for c = C
    g = fullfile(gdir, c.name, 'result.gdf');
    if exist(fullfile(res_dir, [c.name '.mat']), 'file') && exist(g, 'file')
        delete(g); fprintf('  deleted %s (case already processed) to free disk space\n', g);
        f = java.io.File(gdir); if f.getFreeSpace/1e9 >= 10; return; end
    end
end
end
