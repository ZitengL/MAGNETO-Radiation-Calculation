%% make_fig4_bunch_vs_z
% Paper figures: bunch length sigma_z (left axis) and the form factor |F|^2 at lambda = 70 nm of the
% 2000-electron sample (right axis, linear) along the undulator, for
%   fig4_bunch_vs_z_sd0p3   sigma_delta = 0.3 %
%   fig4_bunch_vs_z_en1     epsilon_n = 1 um
% These milder cases replace sigma_delta = 3 % and epsilon_n = 12 um (old figures in figures_unused/), for which
% the bunch is already too long at the undulator entrance and |F|^2 stays at the 1/N = 5e-4 noise level.
% z = mean longitudinal position of the bunch at each GPT output time; range = extent of the undulator field
% map (z = 1 m +- 0.32004 m).
% Data (read only): data/spread_emittance/bunch_vs_z_<case>.mat (compute_bunch_vs_z.m)
% Run: matlab -batch "make_fig4_bunch_vs_z"

this_dir = fileparts(mfilename('fullpath'));
d = fullfile(this_dir, 'data', 'spread_emittance');
zr = 1 + 0.32004*[-1 1];
fs = 8; fs_tick = 7; fname = 'Times New Roman';
c_sz = [0 0 0]; c_F = [0.00 0.45 0.74];
cases = {'sd0p3_en0', '\sigma_\delta = 0.3%', 'fig4_bunch_vs_z_sd0p3'; 'sd0_en1', '\epsilon_n = 1 \mum', 'fig4_bunch_vs_z_en1'};

for c = 1:size(cases, 1)
    S = load(fullfile(d, ['bunch_vs_z_' cases{c, 1} '.mat'])); M = S.M;
    in = M(:, 2) >= zr(1) & M(:, 2) <= zr(2);
    fprintf('%s: sigma_z %.1f -> %.1f nm, |F|^2 %.3g -> %.3g over the undulator\n', cases{c, 1}, ...
        1e9*M(find(in, 1), 3), 1e9*M(find(in, 1, 'last'), 3), M(find(in, 1), 4), M(find(in, 1, 'last'), 4));

    fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 6.0], 'Color', 'w');
    tl  = tiledlayout(fig, 1, 1, 'Padding', 'compact'); ax = nexttile(tl);
    yyaxis(ax, 'left');
    h1 = plot(ax, M(in, 2), 1e9*M(in, 3), '-', 'Color', c_sz, 'LineWidth', 0.9);
    ylabel(ax, '\sigma_z [nm]', 'FontSize', fs); ax.YColor = c_sz; ylim(ax, [0 1.15*max(1e9*M(in, 3))]);
    yyaxis(ax, 'right');
    h2 = plot(ax, M(in, 2), M(in, 4), '-', 'Color', c_F, 'LineWidth', 0.9);
    ylim(ax, [0 0.5]); ax.YColor = c_F;
    ylabel(ax, '|F|^2 at \lambda = 70 nm', 'FontSize', fs);
    set(ax, 'FontName', fname, 'FontSize', fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', 'Box', 'on', ...
        'XLim', zr, 'TickLength', [0.015 0.015]);
    xlabel(ax, '\langlez\rangle [m]', 'FontSize', fs);
    lg = legend(ax, [h1 h2], {'\sigma_z', '|F|^2'}, 'Location', 'northoutside', ...
        'Box', 'off', 'FontName', fname, 'FontSize', fs_tick, 'NumColumns', 2);
    lg.ItemTokenSize = [14 6];
    text(ax, 0.03, 0.80, cases{c, 2}, 'Units', 'normalized', 'FontName', fname, 'FontSize', fs);
    base = fullfile(this_dir, 'figures', cases{c, 3});
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
    exportgraphics(fig, [base '.png'], 'Resolution', 600);
    savefig(fig, [base '.fig']); close(fig);
end
fprintf('Figures written\n');
