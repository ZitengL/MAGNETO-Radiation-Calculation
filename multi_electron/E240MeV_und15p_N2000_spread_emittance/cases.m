function C = cases()
%% cases
% Parameter list of the energy-spread / emittance study (2000 electrons, 15-period undulator, 240 MeV).
%   sig_d : rms relative energy spread sigma_gamma/gamma (uncorrelated, Gaussian)
%   en    : normalised rms emittance, the same in x and y [m rad]; waist at the undulator centre z = 1 m
%           with sigma_x = sigma_y = 20 um there, rms divergence sigma_theta = en/(gamma sigma_x)
% Case names are the folder names.
sd = [0 0.0003 0.001 0.003 0.01 0.03 0     0     0     0     0     0.01];
en = [0 0      0     0     0    0    0.3e-6 1e-6 3e-6 10e-6 12e-6 3e-6];
C = struct('name', {}, 'sig_d', {}, 'en', {});
for i = 1:numel(sd)
    C(i).sig_d = sd(i); C(i).en = en(i);
    C(i).name = sprintf('sd%s_en%s', strrep(sprintf('%g', 100*sd(i)), '.', 'p'), strrep(sprintf('%g', 1e6*en(i)), '.', 'p'));
end
end
