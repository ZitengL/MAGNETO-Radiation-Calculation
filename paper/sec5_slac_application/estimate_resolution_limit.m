%% estimate_resolution_limit
% Sec. 5.3 discussion: up to which sigma_z the calculation resolves the bunch length, as a function of the number of
% macro-particles, and the limit of the setup itself (real 100 pC beam).
% Criterion: the coherent energy in the band equals the noise floor, i.e. the raw energy is twice the floor.
% The floor of N_m macro-particles is nmacro(N_m) times the real incoherent level, nmacro = 6.24e8 / N_m.
% The coherent energy beyond the resolved range is taken from the expected value over the longitudinal positions of
% the same 5e4 macro-particles (E7GeV_und_dipole/results/large_sigz_summary.mat, Wcoh on 10:10:1000 nm), which agrees
% with the direct runs (direct_sigz/check_vs_expected.txt) within 0.7 sd at all eight sigma_z.
% Output: data/resolution_limit.txt. Run: matlab -batch "estimate_resolution_limit"

this_dir = fileparts(mfilename('fullpath'));
src = fullfile(this_dir, '..', '..', 'matlab', 'radiation_multi_electron', 'E7GeV_und_dipole', 'results');
L = load(fullfile(src, 'large_sigz_summary.mat'));
D = load(fullfile(src, 'direct_sigz', 'summary.mat'));
sz = L.sz*1e9; Winc = L.Winc(1:2); nm0 = 12483.02; Ne = 100e-12/1.602176634e-19;
Nm = [5e4 5e5 5e6 5e7 Ne];
lim = zeros(numel(Nm), 2);
for b = 1:2
    lc = log(max(L.Wcoh(b, :), realmin));
    for i = 1:numel(Nm)
        floor_i = Winc(b)*Ne/Nm(i);                      % nmacro * real incoherent
        k = find(L.Wcoh(b, :) < floor_i, 1);
        lim(i, b) = interp1(lc(k-1:k), sz(k-1:k), log(floor_i));
    end
end
fid = fopen(fullfile(this_dir, 'data', 'resolution_limit.txt'), 'w');
fprintf(fid, '%% sigma_z [nm] where the coherent energy equals the macro-particle noise floor (expected value of the 5e4-particle beam)\n');
fprintf(fid, '%% N_macro  nmacro  D1(50-280nm)  D2(300-1000nm); last row = real beam (floor = incoherent level)\n');
for i = 1:numel(Nm); fprintf(fid, '%10.3g %10.3g %8.1f %8.1f\n', Nm(i), Ne/Nm(i), lim(i, :)); end
fprintf(fid, '%% real incoherent level: D1 %.3e, D2 %.3e J/sr; electrons per bunch %.4g; fields of summary.mat: %s\n', Winc, Ne, strjoin(fieldnames(D)', ' '));
fclose(fid);
type(fullfile(this_dir, 'data', 'resolution_limit.txt'));
