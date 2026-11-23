%% compute_polarizer_bands
% Numbers for Sec. 5.1 at sigma_z = 10 nm, (0, -2 mm, 5.9 m), 100 pC: band energies [J/sr] of the undulator alone,
% the dipole alone and both together, without polarizer (x + y + z) and behind the x polarizer (x component),
% in D1 (50-280 nm), D2 (300-1000 nm) and 50-1000 nm, plus the dipole share of the x-polarized energy.
% Two versions, both from E7GeV_und_dipole/results/large_sigz.mat (GPT run at sigma_z = 10 nm, 5e4 macro-particles):
%   EV = expected value over the longitudinal positions, real incoherent level (same as fig_slac_polarizer)
%   DS = direct coherent sum over the macro-particles at their actual positions (R.Er at sz_real = 10 nm,
%        identical to the GPT positions because sz_gpt = 10 nm); contains the macro-particle self term.
% Dipole share in x = (W_x(und + dip) - W_x(und)) / W_x(und + dip), i.e. dipole alone + interference.
% Uncertainty: jackknife over the 25 chunks (100 pC kept: fields x f, self terms x f^2, nmacro x f).
% Output: data/polarizer_bands.txt, data/polarizer_bands.mat. Run: matlab -batch "compute_polarizer_bands"

this_dir = fileparts(mfilename('fullpath'));
src = fullfile(this_dir, '..', '..', 'matlab', 'radiation_multi_electron', 'E7GeV_und_dipole', 'results');
data_dir = fullfile(this_dir, 'data');
c0 = 299792458; eps0 = 8.85418781e-12; sz = 10e-9;
R = load(fullfile(src, 'large_sigz.mat')); R = R.R;
dn = find(R.done); nch = numel(dn); assert(nch == numel(R.done), 'not all chunks done');
assert(abs(R.sz_real(1) - sz) < 1e-12 && abs(R.sz_gpt - sz) < 1e-12, 'sz_real(1) or sz_gpt is not 10 nm');
lam = 1e9*R.lambda(:); om = R.omega(:); k = om/c0; pref = eps0*c0*norm(R.obs)^2/pi; nm = R.nmacro;
F2 = exp(-(k*sz).^2);
bands = [50 280; 300 1000; 50 1000]; bname = {'D1 50-280 nm', 'D2 300-1000 nm', '50-1000 nm'};
Wb = @(S, b) -trapz(om(lam >= bands(b,1) & lam <= bands(b,2)), S(lam >= bands(b,1) & lam <= bands(b,2)));

% per chunk, as (chunk, omega): E0 and Er for part p (1 und, 2 dip, 3 = 1 + 2), component c; self terms I (part 3 = |Eu+Ed|^2)
E0p = @(p, c) squeeze(R.E0(dn, p, c, :));
Erp = @(p, c) squeeze(R.Er(dn, 1, p, c, :));
E0 = cell(3, 3); Er = cell(3, 3); I = cell(3, 3);
for c = 1:3
    E0{1,c} = E0p(1, c); E0{2,c} = E0p(2, c); E0{3,c} = E0{1,c} + E0{2,c};
    Er{1,c} = Erp(1, c); Er{2,c} = Erp(2, c); Er{3,c} = Er{1,c} + Er{2,c};
    for p = 1:3; I{p,c} = squeeze(R.I(dn, p, c, :)); end
end

% quantities: columns 1-3 und, dip, tot (x+y+z); 4-6 und, dip, tot (x only)
nq = 6; qp = [1 2 3 1 2 3]; qc = {1:3, 1:3, 1:3, 1, 1, 1};
Sev = @(p, cs, m, f) sum(cell2mat(arrayfun(@(c) pref*(F2.*(abs(f*sum(E0{p,c}(m,:), 1).').^2 - f^2*sum(I{p,c}(m,:), 1).') ...
                     + f^2*sum(I{p,c}(m,:), 1).'/(nm*f)), cs, 'UniformOutput', false)), 2);
Sds = @(p, cs, m, f) sum(cell2mat(arrayfun(@(c) pref*abs(f*sum(Er{p,c}(m,:), 1).').^2, cs, 'UniformOutput', false)), 2);
Sself = @(p, cs) sum(cell2mat(arrayfun(@(c) pref*sum(I{p,c}, 1).', cs, 'UniformOutput', false)), 2);   % macro-particle self term

meth = {Sev, Sds}; mname = {'EV (expected value, real incoherent level)', 'DS (direct sum at the actual positions, with macro-particle self term)'};
W = zeros(2, 3, nq); Wsd = W; share = zeros(2, 3); share_sd = share;
for im = 1:2
    S = meth{im};
    for b = 1:3
        w = arrayfun(@(q) Wb(S(qp(q), qc{q}, 1:nch, 1), b), 1:nq);
        wj = zeros(nch, nq);
        for c = 1:nch
            m = setdiff(1:nch, c); f = nch/(nch - 1);
            wj(c, :) = arrayfun(@(q) Wb(S(qp(q), qc{q}, m, f), b), 1:nq);
        end
        jsd = @(A) sqrt((nch - 1)/nch*sum((A - mean(A, 1)).^2, 1));
        W(im, b, :) = w; Wsd(im, b, :) = jsd(wj);
        share(im, b) = (w(6) - w(4))/w(6); share_sd(im, b) = jsd((wj(:,6) - wj(:,4))./wj(:,6));
    end
end
Wself = zeros(3, 2);   % macro-particle self term, x component: total, dipole
for b = 1:3; Wself(b, :) = [Wb(Sself(3, 1), b), Wb(Sself(2, 1), b)]; end

fid = fopen(fullfile(data_dir, 'polarizer_bands.txt'), 'w');
fprintf(fid, '%% sigma_z = 10 nm, (0, -2 mm, 5.9 m), 100 pC, 5e4 macro-particles (large_sigz.mat). Energies in J/sr, +- jackknife sd (25 chunks).\n');
fprintf(fid, '%% dipole share in x = (W_x(und+dip) - W_x(und)) / W_x(und+dip)\n');
for im = 1:2
    fprintf(fid, '%% ---- %s\n', mname{im});
    for b = 1:3
        w = squeeze(W(im, b, :)); s = squeeze(Wsd(im, b, :));
        fprintf(fid, '%s\n', bname{b});
        fprintf(fid, '  no polarizer: und %.4g +- %.2g, dip %.4g +- %.2g, und+dip %.4g +- %.2g\n', w(1), s(1), w(2), s(2), w(3), s(3));
        fprintf(fid, '  x polarizer : und %.4g +- %.2g, dip %.4g +- %.2g, und+dip %.4g +- %.2g\n', w(4), s(4), w(5), s(5), w(6), s(6));
        fprintf(fid, '  x(und+dip)/und(no pol) = %.4f, x(und+dip)/(und+dip)(no pol) = %.4f, x(und)/und(no pol) = %.4f\n', w(6)/w(1), w(6)/w(3), w(4)/w(1));
        fprintf(fid, '  dipole share in x = %+.4f +- %.4f\n', share(im, b), share_sd(im, b));
    end
end
fprintf(fid, '%% macro-particle self term, x component (DS only): ');
for b = 1:3; fprintf(fid, '%s total %.3g, dipole %.3g; ', bname{b}, Wself(b, 1), Wself(b, 2)); end
fprintf(fid, '\n');
fclose(fid);
save(fullfile(data_dir, 'polarizer_bands.mat'), 'W', 'Wsd', 'share', 'share_sd', 'Wself', 'bands', 'bname', 'mname');
type(fullfile(data_dir, 'polarizer_bands.txt'));
