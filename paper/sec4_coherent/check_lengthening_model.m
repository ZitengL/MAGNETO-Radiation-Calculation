%% check_lengthening_model
% Ballistic estimate of the coherence loss in the energy-spread / emittance study (on axis).
% The form-factor model with the initial distribution predicts no change. Here each electron gets
% the extra delay it accumulates on the way to the undulator centre (z = 1 m), where its emission
% is represented:
%   energy:  ds_j = -delta_j * L_eff / gamma^2,   L_eff = z_entrance + (1 + K^2/2) * (z_c - z_entrance)
%            (drift from the origin to the undulator entrance, then the slower mean longitudinal
%            velocity inside the undulator up to its centre)
%   angle:   ds_j = +(theta_xj^2 + theta_yj^2) * z_c / 2   (longer path of an electron at an angle)
% and R(omega) = |sum_j exp(-i omega dt_j)|^2 / N with dt_j = (-z_j + ds_j)/(beta c) + observation
% geometry (as in collect_sec4_spread_emittance.m). r_model is weighted with S_inc over the same
% natural-line window as r. This is an estimate (emission is spread over the undulator), not a model
% used in the paper figures.
% Output: data/spread_emittance/lengthening_model.txt
% Run: matlab -batch "check_lengthening_model"

this_dir = fileparts(mfilename('fullpath'));
root = fullfile(this_dir, '..', '..');
gdir = fullfile(root, 'gpt_sims', 'multi_electron', 'E240MeV_und15p_N2000_spread_emittance');
A = load(fullfile(this_dir, 'data', 'spread_emittance', 'sec4_spread_emittance.mat'));
c0 = 299792458; gam0 = 470.66442480; beta0 = sqrt(1 - 1/gam0^2); K = 0.980309;
zc = 1.0; Lu = 15*30.4728e-3; zent = zc - Lu/2; L_eff = zent + (1 + K^2/2)*(zc - zent);
q = 1; on = A.oname{q}; O = A.obs{q}(:)'; hl = A.harm{q}; Wq = A.W.(on); Tq = A.T.(on);
lam = A.lam(:); omega = A.omega(:); names = {A.cases_info.name};
fid = fopen(fullfile(this_dir, 'data', 'spread_emittance', 'lengthening_model.txt'), 'w');
fprintf(fid, '%% On axis. Ballistic lengthening estimate: extra delay -delta*L_eff/gamma^2 (L_eff = %.4f m) + theta^2*%.1f m/2, at the undulator centre.\n', L_eff, zc);
fprintf(fid, '%% r = MAGNETO int S_coh/int S_inc (natural-line window); r_FF = form factor of the initial distribution; r_model = form factor with the extra delays;\n');
fprintf(fid, '%% sig_extra = rms of the extra delay [nm]. Columns: case sig_d en_um n r r_FF r_model sig_extra_nm\n');
fprintf('L_eff = %.4f m\n', L_eff);
for i = 1:numel(names)
    D = readmatrix(fullfile(gdir, names{i}, 'initial_distribution.txt'), 'NumHeaderLines', 1);
    xp = D(:,4)./D(:,6); yp = D(:,5)./D(:,6); dl = D(:,7)/gam0 - 1;
    xw = D(:,1) + xp*zc; yw = D(:,2) + yp*zc; N = numel(dl);
    Lj = vecnorm(O - [xw, yw, zc*ones(N,1)], 2, 2); L0 = norm(O - [0 0 zc]);
    geo = (Lj.^2 - L0^2)./(Lj + L0);
    ds = -dl*L_eff/gam0^2 + (xp.^2 + yp.^2)*zc/2;
    dt = (-D(:,3) + ds)/(beta0*c0) + geo/c0;
    R = zeros(size(omega));
    for k = 1:numel(omega); R(k) = abs(sum(exp(-1i*omega(k)*dt)))^2/N; end
    S = Wq(i).S_inc(:);
    for n = 1:numel(hl)
        w = lam >= Tq.win(i,n,1) & lam <= Tq.win(i,n,2);
        rm = abs(trapz(omega(w), S(w).*R(w)))/abs(trapz(omega(w), S(w)));
        fprintf(fid, '%-12s %-7.4g %-5.3g %d %10.4g %10.4g %10.4g %8.3f\n', names{i}, Wq(i).sig_d, 1e6*Wq(i).en, hl(n), Tq.r_int(i,n), Tq.r_FF_d(i,n), rm, 1e9*std(ds));
        fprintf('%-12s n=%d  r = %9.4g   r_FF = %9.4g   r_model = %9.4g   (sig_extra = %.2f nm)\n', names{i}, hl(n), Tq.r_int(i,n), Tq.r_FF_d(i,n), rm, 1e9*std(ds));
    end
end
fclose(fid);
