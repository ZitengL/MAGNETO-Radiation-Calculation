%% make_distributions
% Initial distributions of all cases in cases.m (one sub-folder per case, with its own GPT input and run.bat).
%   N = 2000 electrons, each macro-particle = 1 electron. All start at t = 0, bunch centre at the origin.
%   x_w, y_w (at the waist z_w = 1 m) and z: the SAME numbers as the baseline run
%     gpt_sims/multi_electron/E240MeV_und15p_N2000 (make_distribution.m: Gaussian, sigma 20 um, 20 um, 10 nm,
%     cut at 4 sigma by redrawing, rng twister seed 20261022)
%   x', y' (divergence) and delta = dgamma/gamma: unit Gaussian deviates from a second stream
%     (seed 20261026, cut at 4 sigma), scaled per case: x' = sigma_theta u1, y' = sigma_theta u2,
%     gamma = gamma0 (1 + sig_d u3), sigma_theta = en/(gamma0 sigma_x). The same deviates in all cases,
%     so that the cases differ only by the parameters.
%   Straight-line drift from the start to the waist: x(0) = x_w - x' z_w, y(0) = y_w - y' z_w.
%   Momentum: GB = sqrt(gamma^2 - 1), GBz = GB/sqrt(1 + x'^2 + y'^2), GBx = x' GBz, GBy = y' GBz.
% Output per case: <case>/initial_distribution.txt and .gdf, <case>/multi_electron.in, <case>/run.bat

clear;
this_dir = fileparts(mfilename('fullpath'));
N = 2000; sig0 = [20e-6 20e-6 10e-9]; ncut = 4; zw = 1.0;
E0 = 240e6; me_g = 9.10953e-31; qe_g = 1.6021892e-19; c0 = 299792458;
G0 = 1 + E0*qe_g/(me_g*c0^2);

% baseline positions (identical to E240MeV_und15p_N2000/make_distribution.m)
rng(20261022, 'twister'); X = zeros(N, 3);
for k = 1:3; X(:, k) = cut_randn(N, ncut) * sig0(k); end
B = readmatrix(fullfile(fileparts(this_dir), 'E240MeV_und15p_N2000', 'initial_distribution.txt'), 'NumHeaderLines', 1);
fprintf('positions identical to the baseline run: max |diff| %.1e m\n', max(abs(X - B(:,1:3)), [], 'all'));
% divergence and energy deviates
rng(20261026, 'twister'); U = zeros(N, 3);
for k = 1:3; U(:, k) = cut_randn(N, ncut); end
fprintf('deviates u1 u2 u3: mean %s, std %s\n', mat2str(mean(U), 3), mat2str(std(U), 4));

tmpl = fileread(fullfile(this_dir, 'multi_electron_template.in'));
bat  = fileread(fullfile(this_dir, 'run_template.bat'));
Cs = cases();
for i = 1:numel(Cs)
    c = Cs(i); d = fullfile(this_dir, c.name); if ~exist(d, 'dir'); mkdir(d); end
    sth = c.en/(G0*sig0(1));
    xp = sth*U(:,1); yp = sth*U(:,2); G = G0*(1 + c.sig_d*U(:,3));
    x0 = X(:,1) - xp*zw; y0 = X(:,2) - yp*zw;
    GB = sqrt(G.^2 - 1); GBz = GB ./ sqrt(1 + xp.^2 + yp.^2);
    T = [x0, y0, X(:,3), xp.*GBz, yp.*GBz, GBz, G, ones(N, 1)];
    fid = fopen(fullfile(d, 'initial_distribution.txt'), 'wt');
    fprintf(fid, 'x y z GBx GBy GBz G nmacro\n'); fprintf(fid, '%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n', T'); fclose(fid);
    [st, msg] = system(sprintf('cd /d "%s" && "%%GPTBIN%%\\asci2gdf.exe" -o initial_distribution.gdf initial_distribution.txt', d));
    if st ~= 0; error('asci2gdf failed: %s', msg); end
    write_crlf(fullfile(d, 'multi_electron.in'), strrep(strrep(tmpl, 'SIGD', sprintf('%g', c.sig_d)), 'ENRM', sprintf('%g', c.en)));
    write_crlf(fullfile(d, 'run.bat'), bat);
    % hard bound of the offset inside the undulator (|z - 1 m| <= 0.32 m), plus the wiggle amplitude ~10 um
    xmax = max(abs(X(:,1:2)) + 0.32*abs([xp yp]), [], 'all') + 10e-6;
    fprintf('%-14s sig_d %-7g en %-6g um: sigma_theta %.4e rad, sample std(delta) %.4e, max |x|,|y| in the undulator <= %.3f mm, sigma_x at the start %.1f um\n', ...
        c.name, c.sig_d, 1e6*c.en, sth, std(G/G0 - 1), 1e3*xmax, 1e6*std(x0));
end

function v = cut_randn(N, ncut)
v = randn(N, 1); out = abs(v) > ncut;
while any(out); v(out) = randn(nnz(out), 1); out = abs(v) > ncut; end
end

function write_crlf(f, t)
t = regexprep(t, '\r?\n', '\r\n');
if any(t > 127); error('non-ASCII character in %s', f); end
fid = fopen(f, 'w'); fwrite(fid, t, 'char'); fclose(fid);
end
