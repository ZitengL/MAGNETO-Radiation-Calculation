%% large_sigz_estimate
% User question (2026-11-05): SLAC says sigma_z may be 10-1000 nm. How does the chosen point (x = 0, y = -2.0 mm,
% z = 5.9 m, x polariser) behave over that range, and how sensitive is it?
% QUICK ESTIMATE (model, not the 5e4 simulation):
%   S1(lambda) = single ideal electron (undulator turned over + dipole, x polarisation), MAGNETO, dz 60 um
%   W(sz) = int_band S1 [N + N(N-1) exp(-k^2 sz^2) exp(-k^2 sr^2 theta^2)] d omega,  N = 100 pC / e, sr = 20 um,
%   theta from the undulator centre (positions only; at 10 nm this model is 30-45 % high vs the 5e4 run, at 50 nm ~ equal)
% Bands: 50-300 nm (current), 0.3-1 um, 1-3 um, 3-10 um. Wavelengths above ~5 um: MAGNETO keeps only the
% acceleration (radiation) term of the Lienard-Wiechert field; not checked there.
% Output: results/large_sigz_estimate.txt, figures/large_sigz_estimate.png/.fig

clear; close all;
this_dir = fileparts(mfilename('fullpath')); root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
c0 = 299792458; eps0 = 8.85418781e-12; N = 100e-12/1.602176634e-19; sr = 20e-6;
obs = [0 -2e-3 5.9]; pref = eps0*c0*norm(obs)^2/pi; th = norm(obs(1:2))/(obs(3) - 0.4);
lam = logspace(log10(50e-9), log10(10e-6), 1500); om = 2*pi*c0./lam; k = om/c0;
tr = read_gpt_trajectory(fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole', 'ideal_electron', 'result_sm_dz60.gdf'));
s = emission_time_spectrum(tr, obs, lam, 'fB', norm(obs)/c0);
S1 = pref*abs(s.Ew(1,:)).^2;                                       % x polarisation [J/(rad/s)/sr] per electron
Ft = exp(-k.^2*sr^2*th^2);
sz = logspace(log10(10e-9), log10(1000e-9), 201);
% bands 2 and 3: SLAC's two detectors (user, 2026-11-05: 50-280 nm and 300-1000 nm)
bands = [50e-9 300e-9; 50e-9 280e-9; 300e-9 1e-6; 1e-6 3e-6; 3e-6 10e-6];
bn = {'50-300 nm', '50-280 nm', '300-1000 nm', '1-3 um', '3-10 um'};
nb = size(bands, 1); Wc = zeros(nb, numel(sz)); Wi = zeros(nb, 1);
for b = 1:nb
    m = lam >= bands(b,1) & lam <= bands(b,2);
    Wi(b) = -trapz(om(m), N*S1(m));
    for i = 1:numel(sz)
        Wc(b,i) = -trapz(om(m), N*(N-1)*S1(m).*Ft(m).*exp(-k(m).^2*sz(i)^2));
    end
end
W = Wc + Wi;
sens = zeros(size(W));
for b = 1:nb; sens(b,:) = -gradient(log(W(b,:)), sz)*1e-9; end      % -dlnW/dsz per nm
fid = fopen(fullfile(this_dir, 'results', 'large_sigz_estimate.txt'), 'w');
fprintf(fid, '%% y = -2.0 mm, x pol., model N S1 + N(N-1) S1 |F|^2 (positions only), W [J/sr]; sens = -dlnW/dsigma_z [1/nm]\n');
fprintf(fid, '%% incoherent level per band: %s\n', strjoin(arrayfun(@(b) sprintf('%s %.2e', bn{b}, Wi(b)), 1:nb, 'UniformOutput', false), ', '));
fprintf(fid, '%% sigma_z_nm'); for b = 1:nb; fprintf(fid, ' | W(%s) coh/incoh sens', bn{b}); end; fprintf(fid, '\n');
for v = [10 20 30 50 70 100 150 200 300 400 500 700 1000]
    [~, i] = min(abs(sz - v*1e-9)); fprintf(fid, '%6.0f', sz(i)*1e9);
    for b = 1:nb; fprintf(fid, ' | %.3e %.2e %.3f', W(b,i), Wc(b,i)/Wi(b), sens(b,i)); end
    fprintf(fid, '\n');
end
for b = 1:nb
    for r = [10 1 0.1]
        i = find(Wc(b,:) < r*Wi(b), 1);
        if isempty(i); fprintf(fid, '%% %s: coherent > %g x incoherent up to 1000 nm\n', bn{b}, r);
        else; fprintf(fid, '%% %s: coherent = %g x incoherent at sigma_z = %.0f nm\n', bn{b}, r, sz(i)*1e9); end
    end
end
% single wavelength (spectrally resolved detector): largest sigma_z with N exp(-k^2 sz^2) |F_perp|^2 = 1
for l = [280 300 1000]*1e-9
    kk = 2*pi/l; fprintf(fid, '%% single lambda %.0f nm: coherent = incoherent at sigma_z = %.0f nm\n', l*1e9, ...
        1e9*sqrt((log(N) - kk^2*sr^2*th^2))/kk);
end
fclose(fid); type(fullfile(this_dir, 'results', 'large_sigz_estimate.txt'));

fh = figure('Position', [100 100 1400 450], 'Visible', 'off');
subplot(1,3,1); loglog(lam*1e9, S1, 'k'); grid on; xlabel('\lambda [nm]'); ylabel('S_1 [J/(rad/s)/sr]');
title('single electron, x pol., y = -2.0 mm');
subplot(1,3,2); for b = 1:nb; loglog(sz*1e9, W(b,:), 'DisplayName', bn{b}); hold on; end
for b = 1:nb; loglog(sz([1 end])*1e9, Wi(b)*[1 1], ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off'); end
grid on; xlabel('\sigma_z [nm]'); ylabel('W in band [J/sr]'); legend('Location', 'southwest'); title('dotted: incoherent level');
subplot(1,3,3); for b = 1:nb; semilogx(sz*1e9, 100*sens(b,:), 'DisplayName', bn{b}); hold on; end
grid on; xlabel('\sigma_z [nm]'); ylabel('-d ln W / d\sigma_z [% per nm]'); legend('Location', 'northwest');
sgtitle('Model estimate (positions only, N = 6.2e8): x = 0, y = -2.0 mm, z = 5.9 m, x polariser');
savefig(fh, fullfile(this_dir, 'figures', 'large_sigz_estimate.fig'));
exportgraphics(fh, fullfile(this_dir, 'figures', 'large_sigz_estimate.png'), 'Resolution', 150);
