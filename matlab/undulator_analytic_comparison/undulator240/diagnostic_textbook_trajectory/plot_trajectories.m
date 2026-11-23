%% plot_trajectories
% x(z) of the electron (240 MeV, single electron) in the three undulator models:
%   (1) real Radia map (paraxial version), GPT: gpt_sims/single_electron/E240MeV_und15p_paraxial
%   (2) ideal field, GPT unduplan (Nu = 16, 1/4, 3/4 end half periods): gpt_sims/single_electron/E240MeV_unduplan_ideal
%   (3) textbook trajectory built in MATLAB (diag_textbook.m): hard-edged By = B0 cos(ku z'), 15 periods,
%       gamma beta_x = K sin(ku z'), x = 0 at the entrance; straight before and after the field
% Panels: x(z) over the whole output range, zooms on the entrance and exit ends, beta_x(z), By(z).
% Output: figures/trajectories_x_vs_z.(png|fig)
clear; close all;
this_dir = fileparts(mfilename('fullpath'));
up_dir = fileparts(this_dir); root_dir = fileparts(fileparts(fileparts(up_dir)));
addpath(fullfile(root_dir, 'magneto'));
c0 = 299792458; me_g = 9.10953e-31; qe_g = 1.6021892e-19;
gam = 470.66442480; K = 0.980309; lamu = 0.030472767; Np = 15; zu = 1.0; ku = 2*pi/lamu;
beta = sqrt(1 - 1/gam^2); B0 = K*2*pi*me_g*c0/(qe_g*lamu);

sdir = fullfile(root_dir, 'gpt_sims', 'single_electron');
R = read_gpt_trajectory(fullfile(sdir, 'E240MeV_und15p_paraxial', 'result.gdf'));
I = read_gpt_trajectory(fullfile(sdir, 'E240MeV_unduplan_ideal', 'result.gdf'));

% textbook trajectory (same construction as diag_textbook.m), extended by straight lines
Lh = Np*lamu/2;
zf = linspace(-Lh, Lh, round(2*Lh/0.25e-6) + 1);
bx = K/gam * sin(ku*zf); bz = sqrt(beta^2 - bx.^2);
xf = cumtrapz(zf, bx ./ bz);
k = 1:40:numel(zf);
if k(end) ~= numel(zf); k = [k numel(zf)]; end
% straight sections before/after; the edge points are repeated so that the field jump is drawn vertically
zt = [0.65, zu - Lh, zu + zf(k), zu + Lh, 1.30]; xt = [0, 0, xf(k), xf(end), xf(end)];
bxt = [0, 0, bx(k), 0, 0]; Byt = [0, 0, B0*cos(ku*zf(k)), 0, 0];

M = {R, I}; lab = {'real map (Radia, paraxial), GPT', 'ideal field (unduplan, Nu = 16), GPT', 'textbook (hard-edged cos, N = 15), MATLAB'};
sty = {'r-', 'b--', 'k:'};
zz = {R.r(3,:), I.r(3,:), zt}; xx = {R.r(1,:), I.r(1,:), xt}; bb = {R.b(1,:), I.b(1,:), bxt}; BB = {R.fB(2,:), I.fB(2,:), Byt};
for j = 1:3
    zin = zz{j} > 0.78 & zz{j} < 1.22;
    [~, iu] = unique(zz{j}, 'last');                   % textbook curve has repeated edge points
    fprintf('%-45s: x range in the field %.3f .. %.3f um, mean x (0.78-1.22 m) %+.3f um, x at z = 1.29 m %+.3f um, beta_x there %+.2e\n', lab{j}, ...
        1e6*min(xx{j}), 1e6*max(xx{j}), 1e6*mean(xx{j}(zin)), 1e6*interp1(zz{j}(iu), xx{j}(iu), 1.29), interp1(zz{j}(iu), bb{j}(iu), 1.29));
end

f = figure('Position', [50 50 1200 950]);
tl = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([1 2]);
for j = 1:3; plot(zz{j}, 1e6*xx{j}, sty{j}, 'LineWidth', 1.1); hold on; end
grid on; xlim([0.65 1.30]); ylabel('x [\mum]'); legend(lab, 'Location', 'southwest');
title('240 MeV single electron: horizontal trajectory x(z) in the three undulator models');
zr = {[0.66 0.86], [1.14 1.30]}; ttl = {'entrance end', 'exit end'};
for p = 1:2
    nexttile;
    for j = 1:3; plot(zz{j}, 1e6*xx{j}, sty{j}, 'LineWidth', 1.1); hold on; end
    grid on; xlim(zr{p}); ylabel('x [\mum]'); title(['x(z), ' ttl{p}]);
end
for p = 1:2
    nexttile;
    for j = 1:3; plot(zz{j}, 1e3*bb{j}, sty{j}, 'LineWidth', 1.1); hold on; end
    grid on; xlim(zr{p}); ylabel('\beta_x [mrad]'); title(['\beta_x(z), ' ttl{p}]);
end
for p = 1:2
    nexttile;
    for j = 1:3; plot(zz{j}, BB{j}, sty{j}, 'LineWidth', 1.1); hold on; end
    grid on; xlim(zr{p}); ylabel('B_y [T]'); xlabel('z [m]'); title(['B_y on the trajectory, ' ttl{p}]);
end
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
if ~exist(fullfile(this_dir, 'figures'), 'dir'); mkdir(fullfile(this_dir, 'figures')); end
exportgraphics(f, fullfile(this_dir, 'figures', 'trajectories_x_vs_z.png'), 'Resolution', 200);
savefig(f, fullfile(this_dir, 'figures', 'trajectories_x_vs_z.fig'));
