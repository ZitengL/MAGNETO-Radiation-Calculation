%% tout_check
% Single electron (undulator turned over + dipole bent up), observers x = 0, y = -1 ... -3 mm, z = 5.9 m:
% output step dz = 120 and 60 um compared with 30 um, per part (undulator / dipole) and polarisation (x / y).
% Output: printed (tout_check.log)
clear;
this_dir = fileparts(mfilename('fullpath')); root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole', 'ideal_electron');
c0 = 299792458; lam = logspace(log10(50), log10(300), 300)*1e-9; zsplit = 0.63;
ys = [-1 -1.5 -2 -2.5 -3]*1e-3; dzs = [30 60 120];
for id = 1:3
    tr = read_gpt_trajectory(fullfile(gdir, sprintf('result_sm_dz%d.gdf', dzs(id))));
    seg = {sub(tr, tr.r(3,:) < zsplit), sub(tr, tr.r(3,:) >= zsplit)};
    for j = 1:numel(ys)
        obs = [0 ys(j) 5.9];
        for is = 1:2; E{id, j, is} = emission_time_spectrum(seg{is}, obs, lam, 'fB', norm(obs)/c0).Ew; end %#ok<SAGROW>
    end
end
pn = {'und', 'dip'}; qn = 'xy';
for id = 2:3
    for j = 1:numel(ys)
        s = sprintf('dz %3d vs 30 um, y = %4.1f mm:', dzs(id), 1e3*ys(j));
        for is = 1:2
            for q = 1:2
                A = abs(E{id,j,is}(q,:)).^2; B = abs(E{1,j,is}(q,:)).^2;
                s = [s sprintf('  %s %s max|dS|/maxS %.1e, band %.1e;', pn{is}, qn(q), max(abs(A-B))/max(B), trapz(A)/trapz(B) - 1)]; %#ok<AGROW>
            end
        end
        % complex field of the x-polarised total (what the polariser passes)
        Ea = E{id,j,1}(1,:) + E{id,j,2}(1,:); Eb = E{1,j,1}(1,:) + E{1,j,2}(1,:);
        s = [s sprintf('  total x: max|dE|/max|E| %.1e', max(abs(Ea - Eb))/max(abs(Eb)))]; %#ok<AGROW>
        disp(s);
    end
end

function t = sub(t, m)
k = find(m); kf = k(any(t.fB(:, k) ~= 0, 1)); k = max(k(1), kf(1)-1) : min(k(end), kf(end)+1);
for f = {'t', 'G'};             t.(f{1}) = t.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; t.(f{1}) = t.(f{1})(:, k); end
end
