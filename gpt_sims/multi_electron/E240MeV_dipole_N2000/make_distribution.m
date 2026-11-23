%% make_distribution
% Initial particle distribution for the multi-electron undulator run.
%   N = 2000 electrons, each macro-particle = 1 electron (q = -e, nmacro = 1)
%   Gaussian in x, y, z: sigma_x = sigma_y = 20 um, sigma_z = 10 nm, truncated at +-4 sigma
%   (points outside are redrawn), pseudo-random (MATLAB rng 'twister', seed 20261022),
%   no energy spread, no angular divergence: GBx = GBy = 0, GBz = sqrt(G^2 - 1),
%   G from 240 MeV kinetic energy with GPT's constants (as in the GPT input files).
%   All particles at t = 0; bunch centre (nominal) at the origin.
% Output: initial_distribution.txt (full precision) and initial_distribution.gdf (asci2gdf).

clear;
this_dir = fileparts(mfilename('fullpath'));

N      = 2000;
sig    = [20e-6, 20e-6, 10e-9];       % sigma_x, sigma_y, sigma_z [m]
ncut   = 4;
seed   = 20261022;
E0     = 240e6;                       % kinetic energy [eV]
me_gpt = 9.10953e-31; qe_gpt = 1.6021892e-19; c0 = 299792458;    % GPT built-in constants
G      = 1 + E0*qe_gpt/(me_gpt*c0^2);
GB     = sqrt(G^2 - 1);

rng(seed, 'twister');
X = zeros(N, 3);
for k = 1:3
    v = randn(N, 1);
    out = abs(v) > ncut;
    while any(out)                    % redraw points beyond the cut
        v(out) = randn(nnz(out), 1);
        out = abs(v) > ncut;
    end
    X(:, k) = v * sig(k);
end

fprintf('G = %.10f, GB = %.10f\n', G, GB);
fprintf('sample mean x, y, z = %.3e, %.3e, %.3e m; sample std = %.4e, %.4e, %.4e m\n', mean(X), std(X));

T = [X, zeros(N, 2), GB*ones(N, 1), G*ones(N, 1), ones(N, 1)];
fid = fopen(fullfile(this_dir, 'initial_distribution.txt'), 'wt');
fprintf(fid, 'x y z GBx GBy GBz G nmacro\n');
fprintf(fid, '%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n', T');
fclose(fid);

cmd = sprintf('cd /d "%s" && "%%GPTBIN%%\\asci2gdf.exe" -o initial_distribution.gdf initial_distribution.txt', this_dir);
[st, msg] = system(cmd);
if st ~= 0; error('asci2gdf failed: %s', msg); end
fprintf('%s\nWrote initial_distribution.txt and initial_distribution.gdf (%d particles)\n', strtrim(msg), N);
