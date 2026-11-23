"""Full precision TXT, comparison statistics and scientific figures for six SRW spectra."""
import analyze as a
import calculate as m
import numpy as np,json
from scipy.io import savemat
from scipy.ndimage import maximum_filter1d
O=m.O
BANDS=[(15,30),(30,60),(60,110),(110,150.00001)]

def sampling(lam,s,dlam,ds):
    assert np.array_equal(lam,dlam[::2]);common=ds[::2]-s;mid=ds[1::2];estimate=(s[:-1]+s[1:])/2
    out={'common_nodes_bitwise_equal':bool(np.array_equal(s,ds[::2])),
        'common_nodes_max_change_over_peak':float(max(abs(common))/max(s)),
        'midpoint_max_error_over_peak':float(max(abs(mid-estimate))/max(s)),
        'peak_baseline_nm':float(lam[np.argmax(s)]),'peak_dense_nm':float(dlam[np.argmax(ds)]),'band_integrals':[]}
    for low,high in BANDS:
        ids=np.flatnonzero((lam>=low*(1-1e-12))&(lam<high));i,j=ids[0],ids[-1]
        v=-np.trapezoid(s[i:j+1],2*np.pi*m.C/(lam[i:j+1]*1e-9))
        w=-np.trapezoid(ds[2*i:2*j+1],2*np.pi*m.C/(dlam[2*i:2*j+1]*1e-9))
        out['band_integrals'].append({'band_nm':[low,min(high,150)],'actual_bounds_nm':[float(lam[i]),float(lam[j])],
            'baseline_J_sr':float(v),'dense_J_sr':float(w),'dense_div_baseline_minus1':float(w/v-1)})
    return out

def run():
    cfg=json.loads((O/'manifest.json').read_text());review=json.loads((O/'pilot_review.json').read_text());ref=np.load(O/'reference.npz')
    lam=ref['lambda_nm'];omega=ref['omega'];ct=review['production_ct_step_um'];ds=review['production_radiation_step_m']
    spectra=[];dense=[];control=[];summaries=[];grid=[];components=[];fields=[];band_rows=[];integral_rows=[];trajectory_rows=[]
    convcols=[np.load(O/'lambda_pilot.npy')];convnames=['lambda_nm']
    for i in range(6):
        base=a.collect('baseline',i,ct,ds);coarse=a.collect('baseline',i,ct,2*ds);fine=a.collect('dense',i,ct,ds)
        s=base['energy'];spectra.append(s);control.append(coarse['energy']);dense.append(fine['energy']);components.append(base['components']);fields.append(base['field'])
        xy=ref['matlab_xy'][i];total=ref['matlab_total'][i];z=ref['matlab_z'][i]
        label=f'{i/10:.1f}';values=np.column_stack([lam,omega,s,*base['components'],xy,total,z,s-xy,s-total,
            s/(m.HBAR*omega),coarse['energy'],s-coarse['energy'],ref['matlab_convergence'][2,i]])
        np.savetxt(O/f'spectrum_theta_{label}mrad.txt',values,fmt='%.17e',header='240MeV single electron, 15-period supplied paraxial undulator map, exact observer D=11*(sin(theta),0,cos(theta)).\n'
            f'theta_x={cfg["theta_rad"][i]:.17e} rad; native SRW gamma and B compensated; original binary and radiation normalization retained.\n'
            'Energy units J/(rad/s)/sr; photon column11 photons/(rad/s)/sr. All5000 original wavelength nodes, no fitting/smoothing/shift.\n'
            'Columns: 1_lambda_nm 2_omega_rad_s 3_SRW_native_xy 4_SRW_x 5_SRW_y 6_MATLAB_original_xy '
            '7_MATLAB_original_total_xyz 8_MATLAB_original_z 9_SRW_minus_MATLAB_xy 10_SRW_minus_MATLAB_total '
            '11_SRW_photons_xy 12_SRW_integration_control_coarser_step 13_SRW_minus_integration_control 14_MATLAB_ct15um_total_xyz\n'
            'MATLAB original spectrum uses ct60um; column14 is its saved fine-step control total, with no unavailable longitudinal reconstruction.')
        dlam=fine['lambda_nm'];np.savetxt(O/f'spectrum_dense_theta_{label}mrad.txt',np.column_stack([dlam,*fine['components'],fine['energy']]),fmt='%.17e',
            header='9999 independently recalculated wavelengths; original nodes plus geometric midpoints, no interpolation to generate this spectrum.\n'
            'Columns: lambda_nm SRW_x SRW_y SRW_xy; energy J/(rad/s)/sr.')
        sg=sampling(lam,s,dlam,fine['energy']);grid.append(sg)
        np.savetxt(O/f'wavelength_midpoints_theta_{label}mrad.txt',np.column_stack([dlam[1::2],fine['energy'][1::2],(s[:-1]+s[1:])/2]),fmt='%.17e',
            header='Sampling diagnostic only; formal dense spectrum is a direct SRW calculation.\nColumns: midpoint_lambda_nm direct_SRW_energy linear_in_loglambda_estimate_from_original_nodes; J/(rad/s)/sr.')
        stats={'theta_mrad':i/10,'SRW_peak_lambda_nm':float(lam[np.argmax(s)]),'SRW_peak_energy':float(max(s)),
            'MATLAB_peak_lambda_nm':float(lam[np.argmax(total)]),'MATLAB_peak_energy':float(max(total)),
            'SRW_vs_MATLAB_xy':a.measure(lam,s,xy),'SRW_vs_MATLAB_total':a.measure(lam,s,total),
            'integration_control':a.measure(lam,s,coarse['energy']),
            'MATLAB60_vs15_total':a.measure(lam,ref['matlab_convergence'][0,i],ref['matlab_convergence'][2,i]),
            'MATLAB30_vs15_total':a.measure(lam,ref['matlab_convergence'][1,i],ref['matlab_convergence'][2,i]),
            'max_MATLAB_z_over_peak':float(max(z)/max(total)),
            'compute_seconds':{'baseline':base['seconds'],'control':coarse['seconds'],'dense':fine['seconds']}}
        ratio=abs(s-xy)/xy;worst=int(np.argmax(ratio))
        stats['worst_pointwise_reference']={'lambda_nm':float(lam[worst]),'MATLAB_energy':float(xy[worst]),
            'MATLAB_energy_over_peak':float(xy[worst]/max(xy)),'signed_relative_difference':float((s[worst]-xy[worst])/xy[worst])}
        summaries.append(stats)
        half=round(.01/np.mean(np.diff(np.log(lam))));env=maximum_filter1d(xy,size=2*half+1,mode='nearest')
        np.savetxt(O/f'relative_errors_theta_{label}mrad.txt',np.column_stack([lam,(s-xy)/xy,(s-xy)/env,(s-xy)/max(xy),(s-coarse['energy'])/env]),fmt='%.17e',
            header='Raw signed relative differences; no intensity floor, threshold, smoothing or excluded samples. Weak minima can have large pointwise ratios.\n'
            'Columns: lambda_nm SRW_minus_MATLABxy_div_MATLABxy SRW_minus_MATLABxy_div_local_envelope SRW_minus_MATLABxy_div_global_peak primary_minus_coarse_SRW_div_MATLAB_local_envelope\n'
            'Local envelope=max within +/-1% logarithmic wavelength. These are diagnostic normalizations, not modifications of the saved spectra.')
        for low,high in BANDS:
            use=(lam>=low*(1-1e-12))&(lam<high)
            band_rows.append([i/10,low,min(high,150),max(abs(s[use]-xy[use])),max(abs(s[use]-xy[use]))/max(xy),
                max(abs(s[use]-xy[use])/env[use]),max(abs(s[use]-coarse['energy'][use])/env[use])])
        for r in sg['band_integrals']:integral_rows.append([i/10,*r['band_nm'],*r['actual_bounds_nm'],r['baseline_J_sr'],r['dense_J_sr'],r['dense_div_baseline_minus1']])
        for d in (60e-6,30e-6,15e-6):
            q=a.collect('pilot',i,7.5,d);convcols.append(q['energy']);convnames.append(f'angle{i/10:g}mrad_ct7p5_rad{d*1e6:g}um')
        if i in (0,3,5):
            pf=a.collect('pilot',i,7.5,15e-6)
            for step in (60,30,15):
                q=a.collect('pilot',i,step,15e-6)
                trajectory_rows.extend(np.column_stack([np.full(len(q['energy']),i/10),q['lambda_nm'],np.full(len(q['energy']),step),q['energy'],pf['energy'],q['energy']-pf['energy']]).tolist())
    S=np.array(spectra);D=np.array(dense);Ctrl=np.array(control);Comp=np.array(components);F=np.array(fields)
    for name,q,unit in [('energy',S,'J/(rad/s)/sr'),('photons',S/(m.HBAR*omega),'photons/(rad/s)/sr')]:
        np.savetxt(O/f'spectra_{name}_all_angles.txt',np.column_stack([lam,omega,q.T]),fmt='%.17e',
            header=f'240MeV single electron, native SRW Ex+Ey; units {unit}. All original5000 wavelengths.\nColumns: lambda_nm omega_rad_s theta0mrad theta0p1mrad theta0p2mrad theta0p3mrad theta0p4mrad theta0p5mrad')
    np.savetxt(O/'pilot_integration_checks.txt',np.column_stack(convcols),fmt='%.17e',header='Sparse native SRW checks; all energy columns J/(rad/s)/sr.\nColumns: '+' '.join(convnames))
    np.savetxt(O/'pilot_trajectory_checks.txt',trajectory_rows,fmt='%.17e',header='Sparse trajectory-sampling checks at 0,0.3,0.5mrad, fixed native radiation integration15um. No formal spectrum replaced.\nColumns: theta_mrad lambda_nm trajectory_ct_step_um native_energy_xy reference_native_ct7p5um_energy_xy energy_difference\nEnergy units J/(rad/s)/sr.')
    np.savetxt(O/'band_comparisons.txt',band_rows,fmt='%.17e',header='Columns: angle_mrad lambda_min_nm lambda_max_nm max_abs_SRW_minus_MATLABxy max_abs_diff_over_full_MATLAB_peak max_diff_over_MATLAB_local_envelope max_fullgrid_integration_diff_over_MATLAB_local_envelope\nLocal envelope: +/-1% logarithmic wavelength window, same convention as original MATLAB. No intensity threshold or smoothing of formal data.')
    np.savetxt(O/'wavelength_band_integrals.txt',integral_rows,fmt='%.17e',header='Columns: angle_mrad nominal_band_min_nm nominal_band_max_nm actual_common_low_nm actual_common_high_nm original_grid_integral_J_sr dense_grid_integral_J_sr dense_div_original_minus1\nBoth integrations in omega, exact same band endpoints. Formal curves evaluated directly.')
    np.savetxt(O/'wavelength_grid_summary.txt',[[i/10,int(r['common_nodes_bitwise_equal']),r['common_nodes_max_change_over_peak'],r['midpoint_max_error_over_peak'],r['peak_baseline_nm'],r['peak_dense_nm']] for i,r in enumerate(grid)],fmt='%.17e',
        header='Columns: theta_mrad common_original_nodes_bitwise_equal_1or0 max_common_node_change_over_peak max_midpoint_reconstruction_error_over_peak peak_lambda_original_nm peak_lambda_dense_nm\nCommon-node agreement alone is not a convergence bound; midpoint estimates are diagnostic only, spectra evaluated directly.')
    np.savetxt(O/'benchmark_summary.txt',[[r['theta_mrad'],r['SRW_peak_lambda_nm'],r['SRW_peak_energy'],r['MATLAB_peak_lambda_nm'],r['MATLAB_peak_energy'],
        r['SRW_vs_MATLAB_xy']['max_difference_over_peak'],r['SRW_vs_MATLAB_xy']['max_difference_over_local_envelope'],r['integration_control']['max_difference_over_local_envelope']] for r in summaries],fmt='%.17e',
        header='Columns: angle_mrad SRW_peak_lambda_nm SRW_peak_energy MATLAB_peak_lambda_nm MATLAB_peak_energy max_SRW_minus_MATLABxy_over_MATLABpeak max_SRW_minus_MATLABxy_over_MATLAB_local_envelope max_integration_difference_over_native_local_envelope\nEnergy J/(rad/s)/sr. Raw differences; no wavelength/normalization fit.')
    theta=ref['theta_rad'];np.savetxt(O/'harmonic_markers.txt',np.column_stack([theta*1e3,np.arctan(11*np.sin(theta)/(11*np.cos(theta)-1))*1e3,ref['harmonics_nm']]),fmt='%.17e',
        header='Original MATLAB theoretical markers, not fitted SRW peak locations.\nColumns: theta_origin_mrad theta_seen_from_undulator_center_mrad lambda_n1_nm lambda_n2_nm lambda_n3_nm lambda_n4_nm lambda_n5_nm')
    arrays={'lambda_nm':lam,'omega':omega,'theta_rad':theta,'SRW_energy_xy':S,'SRW_components':Comp,'SRW_field_native':F,
        'SRW_photons_xy':S/(m.HBAR*omega),'MATLAB_xy':ref['matlab_xy'],'MATLAB_total':ref['matlab_total'],'MATLAB_z':ref['matlab_z'],
        'SRW_integration_control':Ctrl,'dense_lambda_nm':dlam,'dense_SRW_xy':D}
    np.savez_compressed(O/'six_spectra_results.npz',**arrays);savemat(O/'six_spectra_results.mat',arrays,do_compression=True)
    summary={'status':'all native SRW spectra complete','angles':summaries,'wavelength_grid':grid,
        'settings':{'ct_step_um':ct,'radiation_step_m':ds,'control_radiation_step_m':2*ds,'original_nodes':5000,'dense_nodes':9999},
        'timing':json.loads((O/'full_completed.json').read_text()),'pilot':review}
    (O/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    from plots import draw
    draw(arrays,ref,summaries)
    verify(arrays)
    print(json.dumps({'angles':summaries,'wall_seconds':summary['timing']['wall_seconds']},indent=2),flush=True)

def verify(arrays):
    for i in range(6):
        q=np.loadtxt(O/f'spectrum_theta_{i/10:.1f}mrad.txt');assert q.shape==(5000,14) and np.all(np.isfinite(q))
        assert np.array_equal(q[:,2],arrays['SRW_energy_xy'][i]);assert max(abs(q[:,2]-q[:,3]-q[:,4]))/max(q[:,2])<1e-14
        assert max(abs(q[:,10]*(m.HBAR*q[:,1])-q[:,2]))/max(q[:,2])<1e-14
        d=np.loadtxt(O/f'spectrum_dense_theta_{i/10:.1f}mrad.txt');assert d.shape==(9999,4) and np.all(np.isfinite(d))
        assert np.array_equal(d[:,3],arrays['dense_SRW_xy'][i])
    E=np.loadtxt(O/'spectra_energy_all_angles.txt');N=np.loadtxt(O/'spectra_photons_all_angles.txt')
    assert E.shape==N.shape==(5000,8) and np.all(np.isfinite(E)) and np.all(np.isfinite(N))
    assert np.array_equal(E[:,2:].T,arrays['SRW_energy_xy'])
    (O/'txt_verified.json').write_text(json.dumps({'native_TXT_roundtrip_exact':True,'finite_data':True,'components_and_photon_conversion':True,'angles':6,'original_nodes':5000,'dense_nodes':9999},indent=2))

if __name__=='__main__':run()
