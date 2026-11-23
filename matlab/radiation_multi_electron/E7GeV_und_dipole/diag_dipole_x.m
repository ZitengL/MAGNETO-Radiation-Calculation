%% diag_dipole_x
% Is the dipole x-polarised light in the 5e4 run (results/combined.mat) physical or macro-particle noise?
%  - W_dip_x (coherent sum of the dipole part) vs the dipole self terms sum_j |E_dip,x,j|^2: if they are equal and
%    W_dip_x does not change with sigma_z, it is the random-phase noise of the macro-particles
%    (the real incoherent level is self / nmacro, nmacro = N/5e4).
%  - interference I = W_tot_x - W_und_x - W_dip_x with its jackknife error over the 25 chunks.
% Output: results/diag_dipole_x.txt

clear;
this_dir = fileparts(mfilename('fullpath')); res_dir = fullfile(this_dir, 'results');
c0 = 299792458; eps0 = 8.85418781e-12;
R = load(fullfile(res_dir, 'combined.mat')); R = R.R;
om = R.omega; szl = R.sigma_z; no = size(R.obs, 1); nch = size(R.Ec, 1);
W = @(S) -trapz(om, S(:)');
Iall = squeeze(sum(R.Ic, 1));                                      % obs x part x comp x omega
fid = fopen(fullfile(res_dir, 'diag_dipole_x.txt'), 'w');
fprintf(fid, '%% band 50-300 nm [J/sr], x polarisation; nmacro = %.1f\n', R.nmacro);
fprintf(fid, '%% y_mm sigma_z_nm W_und_x W_dip_x W_dip_x_self W_und_x_self interference interference_jk_sd contamination contamination_jk_sd\n');
for o = 1:no
    pref = eps0*c0*norm(R.obs(o,:))^2/pi;
    sd = W(pref*squeeze(Iall(o,2,1,:))); su = W(pref*squeeze(Iall(o,1,1,:)));
    for s = 1:numel(szl)
        Eu = squeeze(R.Ec(:,s,o,1,1,:)); Ed = squeeze(R.Ec(:,s,o,2,1,:));   % chunk x omega
        f = @(Eu, Ed) [W(pref*abs(sum(Eu,1)).^2), W(pref*abs(sum(Ed,1)).^2), W(pref*abs(sum(Eu+Ed,1)).^2)];
        v = f(Eu, Ed); I = v(3) - v(1) - v(2); C = (v(3) - v(1))/v(1);
        jI = zeros(1, nch); jC = jI;
        for j = 1:nch
            m = setdiff(1:nch, j); w = f(Eu(m,:)*nch/(nch-1), Ed(m,:)*nch/(nch-1));
            jI(j) = w(3) - w(1) - w(2); jC(j) = (w(3) - w(1))/w(1);
        end
        sI = sqrt((nch-1)/nch*sum((jI - mean(jI)).^2)); sC = sqrt((nch-1)/nch*sum((jC - mean(jC)).^2));
        fprintf(fid, '%5.1f %3.0f %.4e %.4e %.4e %.4e %+.4e %.2e %+.4e %.2e\n', 1e3*R.obs(o,2), 1e9*szl(s), v(1), v(2), sd, su, I, sI, C, sC);
    end
end
fclose(fid); type(fullfile(res_dir, 'diag_dipole_x.txt'));
