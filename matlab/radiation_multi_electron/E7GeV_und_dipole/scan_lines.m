%% scan_lines
% Ideal electron through undulator + dipole (bent up), z = 5.9 m plane, single electron (no beam averaging):
% energy in 50-300 nm of the undulator part, the dipole part and the sum, per polarisation (x, y),
% along the x axis (y = 0, wiggle plane of the undulator) and the y axis (x = 0, bend plane of the dipole).
% Both undulator polarities (sp: as before, sm: turned over) because the real relative polarity is not known.
% Also the same energies weighted with exp(-k^2 sz^2), sz = 10 nm (longitudinal form factor only, a rough
% indicator of which wavelengths count for the coherent signal; NOT the full coherent calculation).
% Output: results/scan_lines.mat, results/scan_lines_{x,y}.txt, figures/scan_lines.png

clear; close all;
this_dir = fileparts(mfilename('fullpath')); root_dir = fileparts(fileparts(fileparts(this_dir)));
addpath(fullfile(root_dir, 'magneto'));
gdir = fullfile(root_dir, 'gpt_sims', 'multi_electron', 'E7GeV_und_dipole', 'ideal_electron');
res_dir = fullfile(this_dir, 'results'); fig_dir = fullfile(this_dir, 'figures');
c0 = 299792458; eps0 = 8.85418781e-12;
lam = logspace(log10(50), log10(300), 300)*1e-9; om = 2*pi*c0./lam; k = om/c0; Fz = exp(-k.^2*(10e-9)^2);
xs = [0.1 0.25:0.25:3 3.5:0.5:6 7:1:20]*1e-3; ys = [-10:1:-4 -3.5:0.25:-0.25 -0.1 0.1 0.25:0.25:3.5 4:1:10]*1e-3;
lines = {[xs(:) zeros(numel(xs),1)], [zeros(numel(ys),1) ys(:)]}; lname = {'x', 'y'};
pol = {'sp', 'sm'}; zsplit = 0.63;
for ip = 1:2
    tr = read_gpt_trajectory(fullfile(gdir, ['result_' pol{ip} '_dz60.gdf']));
    seg = {sub(tr, tr.r(3,:) < zsplit), sub(tr, tr.r(3,:) >= zsplit)};
    for il = 1:2
        P = lines{il};
        for j = 1:size(P, 1)
            obs = [P(j,:) 5.9]; p = eps0*c0*norm(obs)^2/pi;
            Eu = emission_time_spectrum(seg{1}, obs, lam, 'fB', norm(obs)/c0).Ew;
            Ed = emission_time_spectrum(seg{2}, obs, lam, 'fB', norm(obs)/c0).Ew;
            for q = 1:2                                % 1: x pol, 2: y pol
                for w = 1:2                            % 1: plain, 2: weighted with Fz
                    g = p*(w == 1) + p*Fz*(w == 2);
                    S.(lname{il}).Wu(j, q, w, ip) = -trapz(om, g.*abs(Eu(q,:)).^2);
                    S.(lname{il}).Wd(j, q, w, ip) = -trapz(om, g.*abs(Ed(q,:)).^2);
                    S.(lname{il}).Wt(j, q, w, ip) = -trapz(om, g.*abs(Eu(q,:) + Ed(q,:)).^2);
                end
            end
        end
    end
end
S.x.pos = xs; S.y.pos = ys;
save(fullfile(res_dir, 'scan_lines.mat'), 'S', 'lam');
for il = 1:2
    L = S.(lname{il}); n = numel(L.pos);
    fid = fopen(fullfile(res_dir, sprintf('scan_lines_%s.txt', lname{il})), 'wt');
    fprintf(fid, '%% single electron, z = 5.9 m, %s axis, 50-300 nm energy [J/sr], undulator polarity sm (turned over) and sp\n', lname{il});
    fprintf(fid, '%% pos_mm  [sm] Wund_x Wund_y Wdip_x Wdip_y Wtot_x Wtot_y  [sp] Wtot_x Wtot_y  [sm, Fz(10nm)-weighted] Wund_x Wund_y Wdip_x Wdip_y\n');
    for j = 1:n
        fprintf(fid, '%8.3f %s\n', 1e3*L.pos(j), sprintf('%.4e ', L.Wu(j,1,1,2), L.Wu(j,2,1,2), L.Wd(j,1,1,2), L.Wd(j,2,1,2), L.Wt(j,1,1,2), L.Wt(j,2,1,2), ...
            L.Wt(j,1,1,1), L.Wt(j,2,1,1), L.Wu(j,1,2,2), L.Wu(j,2,2,2), L.Wd(j,1,2,2), L.Wd(j,2,2,2)));
    end
    fclose(fid);
end
f = figure('Position', [50 50 1200 800]); cs = {'x', 'y'};
for il = 1:2
    L = S.(lname{il}); u = 1e3*L.pos;
    for q = 1:2
        subplot(2,2,(il-1)*2+q); hold on;
        plot(u, L.Wu(:,q,1,2), 'b.-', 'DisplayName', 'undulator');
        plot(u, L.Wd(:,q,1,2), 'r.-', 'DisplayName', 'dipole');
        plot(u, abs(L.Wt(:,q,1,2) - L.Wu(:,q,1,2)), 'k.-', 'DisplayName', '|total - undulator|, turned over');
        plot(u, abs(L.Wt(:,q,1,1) - L.Wu(:,q,1,1)), 'k.:', 'DisplayName', '|total - undulator|, as before');
        set(gca, 'YScale', 'log'); grid on; xlabel(sprintf('%s at z = 5.9 m [mm]', lname{il}));
        ylabel('single-electron energy, 50-300 nm [J/sr]'); title(sprintf('%s axis, %s polarisation', lname{il}, cs{q}), 'Interpreter', 'none');
        legend('Location', 'best');
    end
end
set(findall(f, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(f, fullfile(fig_dir, 'scan_lines.png'), 'Resolution', 150); close(f);

function t = sub(t, m)
k = find(m); kf = k(any(t.fB(:, k) ~= 0, 1)); k = max(k(1), kf(1)-1) : min(k(end), kf(end)+1);
for f = {'t', 'G'};             t.(f{1}) = t.(f{1})(k);    end
for f = {'r', 'b', 'fE', 'fB'}; t.(f{1}) = t.(f{1})(:, k); end
end
