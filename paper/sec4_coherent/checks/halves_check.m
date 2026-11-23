% Is the two-peak coherent spectrum (sd1, en10) deterministic or a finite-N random pattern?
% Uses the odd/even half sums stored by run_case.m.  E_half = D/2 + R_half (D: deterministic coherent part,
% R: random part of each half, independent).  4 Re(E_o conj(E_e)) estimates |D|^2 (random part averages out);
% |E_o - E_e|^2 contains only random parts and estimates the incoherent sum I_all.
rd = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'matlab', 'radiation_multi_electron', 'und15p_N2000_spread_emittance', 'results');
eps0 = 8.8541878128e-12; c0 = 299792458;
cs = {'sd0_en0','sd0p3_en0','sd1_en0','sd3_en0','sd0_en3','sd0_en10','sd0_en12'};
lam1 = 102.002; out = [fileparts(mfilename('fullpath')) filesep];
fig = figure('Units','centimeters','Position',[1 1 30 22],'Color','w'); tl = tiledlayout(fig, 3, 3);
for c = 1:numel(cs)
    S = load(fullfile(rd, [cs{c} '.mat'])); r = S.res; K = eps0*c0*121/pi;
    x = r.lambda_nm(:)/lam1 - 1; sel = abs(x) < 0.3; win = abs(x) < 0.1;
    Ea = squeeze(r.E(1,1,:,:)); Eo = squeeze(r.E(1,2,:,:)); Ee = squeeze(r.E(1,3,:,:));
    Sa = K*sum(abs(Ea).^2,1)'; So = K*sum(abs(Eo).^2,1)'; Se = K*sum(abs(Ee).^2,1)';
    Sx = K*4*sum(real(Eo.*conj(Ee)),1)'; Sd = K*sum(abs(Eo-Ee).^2,1)'; Ia = K*squeeze(sum(r.I(1,1,:,:),3)); Ia = Ia(:);
    fprintf('%-10s  int over |x|<0.1:  coh %.3g  inc %.3g  4Re(EoEe*) %.3g  |Eo-Ee|^2 %.3g   ratio det/inc %.3g, diff/inc %.3f\n', cs{c}, ...
        trapz(x(win),Sa(win)), trapz(x(win),Ia(win)), trapz(x(win),Sx(win)), trapz(x(win),Sd(win)), ...
        trapz(x(win),Sx(win))/trapz(x(win),Ia(win)), trapz(x(win),Sd(win))/trapz(x(win),Ia(win)));
    ax = nexttile(tl); hold(ax,'on');
    plot(ax, x(sel), Sa(sel), 'k-', x(sel), 4*So(sel), 'r-', x(sel), 4*Se(sel), 'b-', x(sel), Sx(sel), 'g--', x(sel), Ia(sel), 'm:', 'LineWidth', 1);
    title(ax, strrep(cs{c},'_','\_')); xlim(ax, [-0.3 0.3]); grid(ax,'on');
    if c == 1; legend(ax, {'|E_{all}|^2','4|E_{odd}|^2','4|E_{even}|^2','4Re(E_oE_e^*)','\Sigma|E|^2'}, 'Location','northeast'); end
end
exportgraphics(fig, [out 'halves_check.png'], 'Resolution', 110);
