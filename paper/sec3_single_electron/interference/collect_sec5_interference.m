%% collect_sec5_interference
% Collects the undulator + downstream dipole interference results (single
% electron, 240 MeV, 213 distances D) into ./data and copies the original
% figures into ./figures_original. Nothing is recomputed.
%
% Source (read only): matlab/interference_und_dipole/results/interference_scan.mat, D0p5_*.txt
% Run: matlab -batch "collect_sec5_interference"

this_dir = fileparts(mfilename('fullpath'));
root     = fullfile(this_dir, '..', '..', '..');
data_dir = fullfile(this_dir, 'data');             if ~exist(data_dir, 'dir'); mkdir(data_dir); end
fig_dir  = fullfile(this_dir, 'figures_original'); if ~exist(fig_dir, 'dir');  mkdir(fig_dir);  end
src      = fullfile(root, 'matlab', 'interference_und_dipole');

s = load(fullfile(src, 'results', 'interference_scan.mat'));  A = s.A;
lam = A.lambda_nm(:); omega = A.omega(:);
obs = {'onaxis', 'th0p5'};
obs_text = {'on axis (0, 0, 11) m', '11 m*(sin 0.5 mrad, 0, cos 0.5 mrad), horizontal'};

fid_s = fopen(fullfile(data_dir, 'summary_numbers.txt'), 'w');
maps = struct();
for k = 1:2
    o = A.(obs{k});
    D = o.D(:);

    % spectra at D = 0.5 m (copied as written by analyze_scan.m)
    copyfile(fullfile(src, 'results', ['D0p5_' obs{k} '.txt']), fullfile(data_dir, ['spectra_D0p5m_' obs{k} '.txt']));

    % fringe delay vs D (stored in s in the source; written as c*dtau in um)
    cu = 299792458 * 1e6;
    M = [D, cu * o.dtau(:), cu * o.dtau_pred(:), cu * (o.dtau(:) - o.dtau_pred(:))];
    write_txt(fullfile(data_dir, ['delay_vs_D_' obs{k} '.txt']), M, { ...
        sprintf('Undulator + dipole, single electron 240 MeV, observer %s', obs_text{k}), ...
        'dtau = change of the delay of the dipole light relative to the undulator light, measured from the phase of Ex*conj(Ey), relative to D = 0.5 m', ...
        'dtau_pred = geometric prediction from the actual trajectory and the exact observer geometry. Units: c*dtau in um.', ...
        sprintf('Linear fit: measured slope c*dtau/dD = %.7e, predicted %.7e', o.slope_meas, o.slope_pred), ...
        'Columns: 1 D_m  2 c*dtau_meas_um  3 c*dtau_pred_um  4 meas-pred_um'});
    fprintf(fid_s, '%s: slope measured %.7e, predicted %.7e (rel. %.1e); max |meas - pred| over %d D = %.2e um\n', ...
        obs{k}, o.slope_meas, o.slope_pred, o.slope_meas / o.slope_pred - 1, numel(D), cu * max(abs(o.dtau - o.dtau_pred)));

    % delay spectrum |sum Ex Ey* e^{-i w t} dw| for four D
    M = [cu * o.delay_t(:), o.delay_spec'];
    write_txt(fullfile(data_dir, ['delay_spectrum_' obs{k} '.txt']), M, { ...
        sprintf('Delay spectrum |sum_w Ex Ey* exp(-i w t) dw| (arbitrary units, as in the source), observer %s', obs_text{k}), ...
        sprintf('Columns: 1 c*t_um  2-5 D = %s m; box centre c*t = %s um', ...
            strjoin(compose('%.3g', o.delay_D), ', '), strjoin(compose('%.4f', cu * o.delay_peak), ', '))});

    % (lambda, D) maps
    I_rel = o.I_int ./ o.S0;
    maps.(obs{k}) = struct('D_m', D, 'lambda_nm', lam, 'cos_dphi', o.cos_dphi, 'I_int_over_S0', I_rel, ...
        'visibility45', o.visibility45);
    write_map(fullfile(data_dir, ['map_cos_dphi_' obs{k} '.txt']), D, lam, o.cos_dphi, ...
        sprintf('cos(dphi) = S2/(2 sqrt(Sx Sy)), phase difference between Ex and Ey, observer %s', obs_text{k}));
    write_map(fullfile(data_dir, ['map_I_int_over_S0_' obs{k} '.txt']), D, lam, I_rel, ...
        sprintf('(S0 - S_und - S_dip)/S0, interference term in the total intensity, observer %s', obs_text{k}));

    % D = 0.5 m numbers
    [~, i5] = min(abs(D - 0.5));
    fprintf(fid_s, '   D = 0.5 m: median S_dip/S_und = %.2e; max |I_int|/S0 = %.2e (median %.2e); median sqrt(S2^2+S3^2)/S0 = %.2f\n', ...
        median(o.S_dip(i5, :) ./ o.S_und), max(abs(I_rel(i5, :))), median(abs(I_rel(i5, :))), ...
        median(sqrt(o.S2(i5, :).^2 + o.S3(i5, :).^2) ./ o.S0(i5, :)));
    fprintf(fid_s, '   delay-spectrum box centre c*t at D = %s m: %s um\n', ...
        strjoin(compose('%.3g', o.delay_D), ', '), strjoin(compose('%.4f', cu * o.delay_peak), ', '));
end
fclose(fid_s);
save(fullfile(data_dir, 'sec5_interference_maps.mat'), 'maps', '-v7.3');
type(fullfile(data_dir, 'summary_numbers.txt'));

% original figures
f = dir(fullfile(src, 'figures', '*.png'));
for k = 1:numel(f); copyfile(fullfile(src, 'figures', f(k).name), fullfile(fig_dir, f(k).name)); end

function write_txt(fn, M, hdr)
    fid = fopen(fn, 'w');
    for k = 1:numel(hdr); fprintf(fid, '%% %s\n', hdr{k}); end
    fmt = [repmat('%.17e ', 1, size(M, 2) - 1) '%.17e\n'];
    fprintf(fid, fmt, M');
    fclose(fid);
end

function write_map(fn, D, lam, Z, title_line)
    % first data row: NaN then the wavelengths; following rows: D then Z(D, :)
    fid = fopen(fn, 'w');
    fprintf(fid, '%% %s\n', title_line);
    fprintf(fid, '%% Matrix layout: row 1 = [NaN, lambda_nm (5000 values)]; rows 2.. = [D_m, values at each lambda]; 7 significant digits\n');
    fmt = [repmat('%.6e ', 1, numel(lam)) '%.6e\n'];
    fprintf(fid, fmt, [NaN; lam(:)]);
    fprintf(fid, fmt, [D(:), Z]');
    fclose(fid);
end
