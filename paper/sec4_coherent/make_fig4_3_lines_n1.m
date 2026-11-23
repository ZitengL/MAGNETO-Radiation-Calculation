%% make_fig4_3_lines_n1
% Paper figure: n = 1 line shapes on axis with energy spread and emittance, two columns (17 cm), 2 x 2:
%   (a) coherent |sum_j E_j|^2, energy spread      (b) incoherent sum_j |E_j|^2, energy spread
%   (c) coherent, emittance                         (d) incoherent, emittance
% The two columns use different parameters, because coherence is lost long before the incoherent line broadens:
%   coherent   : sigma_delta = 0, 0.1, 0.3 %;  epsilon_n = 0, 1, 3 um   (all still coherent, r >= 58; for larger
%                values the coherent spectrum is a finite-N speckle pattern, see checks/halves_check.m)
%                normalized to the peak of the cold bunch, so the drop of the peak is visible; x zoomed to +-0.15
%   incoherent : sigma_delta = 0, 1, 3 %;  epsilon_n = 0, 10, 12 um, normalized to the incoherent peak of the
%                cold bunch, as the coherent column (user, 2026-10-30)
% Peaks = maximum within |lambda/lambda_1 - 1| < 0.1, lambda_1 = n = 1 peak of the cold incoherent sum.
% Legend: peak value d2W/domega dOmega [J s sr^-1]; for the coherent panels also the effective form factor
%   |F|^2_eff = (r - 1)/(N - 1), r = int S_coh domega / int S_inc domega over the natural-line window (table T.r_int).
% Data (read only): data/spread_emittance/sec4_spread_emittance.mat (collect_sec4_spread_emittance.m)
% Variant (user, 2026-10-30): every curve normalized to its own peak, all four panels ->
%   figures/fig4_3_lines_n1_own_peak.  Run: matlab -batch "norm_mode='own'; make_fig4_3_lines_n1"
% Variant (user, 2026-10-30): normalized to the cold-bunch peak as the main figure, log y axis ->
%   figures/fig4_3_lines_n1_log.  Run: matlab -batch "ylog=true; make_fig4_3_lines_n1"
% Run: matlab -batch "make_fig4_3_lines_n1"

this_dir = fileparts(mfilename('fullpath'));
if ~exist('norm_mode', 'var'); norm_mode = 'cold'; end   % 'cold': cold-bunch peak of each column; 'own': own peak
own = strcmp(norm_mode, 'own');
if ~exist('ylog', 'var'); ylog = false; end
A = load(fullfile(this_dir, 'data', 'spread_emittance', 'sec4_spread_emittance.mat'));
W = A.W.onaxis; T = A.T.onaxis; C = A.cases_info; N = A.N; lam = A.lam(:);
name2idx = @(s) find(strcmp({C.name}, s));
lam1 = W(1).lpk_inc(1);                          % cold incoherent n = 1 peak
x = lam / lam1 - 1;
win = abs(x) < 0.1;
fprintf('lambda_1 = %.3f nm\n', lam1);

fs = 8; fs_tick = 7; fname = 'Times New Roman';
cols = [0 0 0; 0.00 0.45 0.74; 0.85 0.33 0.10];   % cold, middle, largest
% {kind, case names, labels} per panel, in the order (a) (b) (c) (d)
P = {'coh', {'sd0_en0', 'sd0p1_en0', 'sd0p3_en0'}, {'\sigma_\delta = 0', '\sigma_\delta = 0.1%', '\sigma_\delta = 0.3%'}; ...
     'inc', {'sd0_en0', 'sd1_en0',   'sd3_en0'},   {'\sigma_\delta = 0', '\sigma_\delta = 1%',   '\sigma_\delta = 3%'}; ...
     'coh', {'sd0_en0', 'sd0_en1',   'sd0_en3'},   {'\epsilon_n = 0', '\epsilon_n = 1 \mum',  '\epsilon_n = 3 \mum'}; ...
     'inc', {'sd0_en0', 'sd0_en10',  'sd0_en12'},  {'\epsilon_n = 0', '\epsilon_n = 10 \mum', '\epsilon_n = 12 \mum'}};
S0.coh = max(W(name2idx('sd0_en0')).S_coh(win)); % cold-bunch peaks: normalization of each column
S0.inc = max(W(name2idx('sd0_en0')).S_inc(win));

fig = figure('Units', 'centimeters', 'Position', [2 2 17 11], 'Color', 'w', 'Renderer', 'painters');   % no OpenGL (PNG export ran out of memory / broke on 2026-10-30)
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
letters = 'abcd';
for p = 1:4
    kind = P{p, 1}; coh = strcmp(kind, 'coh');
    if coh; xl = 0.15; else; xl = 0.3; end
    sel = abs(x) < xl;
    ax = nexttile(tl); hold(ax, 'on');
    hh = gobjects(1, 3); lab = cell(1, 3);
    for j = 1:3
        k = name2idx(P{p, 2}{j});
        S = W(k).(['S_' kind])(:);
        pk = max(S(win));
        if own; nrm = pk; else; nrm = S0.(kind); end
        hh(j) = plot(ax, x(sel), S(sel) / nrm, '-', 'Color', cols(j, :), 'LineWidth', 0.8);
        if coh
            F2 = (T.r_int(k, 1) - 1) / (N - 1);
            lab{j} = sprintf('%s: %s, |F|^2_{eff} = %.2f', P{p, 3}{j}, sci(pk), F2);
            fprintf('%s coh: peak %.3e (%.3f of cold), r = %.4g, |F|^2_eff = %.3e\n', C(k).name, pk, pk/nrm, T.r_int(k, 1), F2);
        else
            lab{j} = sprintf('%s: %s', P{p, 3}{j}, sci(pk));
            fprintf('%s inc: peak %.3e (%.3f of cold)\n', C(k).name, pk, pk/nrm);
        end
    end
    set(ax, 'FontName', fname, 'FontSize', fs_tick, 'LineWidth', 0.5, 'TickDir', 'in', 'Box', 'on', ...
        'XLim', xl*[-1 1], 'YLim', [0 2.15], 'YTick', 0:0.5:1, 'TickLength', [0.015 0.015]);   % legend above the curves
    if ylog; set(ax, 'YScale', 'log', 'YLim', [1e-2 10^2.6], 'YTick', 10.^(-2:0), 'YMinorTick', 'on'); end   % floor just below the lowest main peak (0.025)
    if coh; ttl = 'coherent |\Sigma_jE_j|^2'; else; ttl = 'incoherent \Sigma_j|E_j|^2'; end
    lg = legend(ax, hh, lab, 'Location', 'north', 'Box', 'off', 'FontName', fname, 'FontSize', fs_tick - 1);
    lg.ItemTokenSize = [12 6];
    lg.Title.String = [ttl ', peak d^2W/d\omegad\Omega [J s sr^{-1}]'];
    lg.Title.FontWeight = 'normal'; lg.Title.FontSize = fs_tick - 1; lg.Title.FontName = fname;
    if p > 2; xlabel(ax, '\lambda/\lambda_1 - 1', 'FontSize', fs); end
    if own; ylabel(ax, 'S / S_{peak}', 'FontSize', fs); else; ylabel(ax, 'S / S_{peak}(cold)', 'FontSize', fs); end
    text(ax, 0.03, 0.95, ['(' letters(p) ')'], 'Units', 'normalized', 'FontName', fname, 'FontSize', fs, 'VerticalAlignment', 'top');
end

base = fullfile(this_dir, 'figures', 'fig4_3_lines_n1');
if own; base = [base '_own_peak']; end
if ylog; base = [base '_log']; end
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [base '.png'], 'Resolution', 600);
savefig(fig, [base '.fig']); close(fig);
fprintf('Written %s\n', base);

function s = sci(v)
% a.b x 10^n in TeX
e = floor(log10(abs(v))); m = v / 10^e;
if round(m, 1) >= 10; m = m / 10; e = e + 1; end
s = sprintf('%.1f\\times10^{%d}', m, e);
end
