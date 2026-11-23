"""Exact TXT exports and static scientific figures; no fitted correction or smoothing."""
import sys,json
sys.dont_write_bytecode=True
import results as r
import model as m
import numpy as np
O=m.O;summary=json.loads((O/'summary.json').read_text());cfg=json.loads((O/'manifest.json').read_text())
verified={};refinement={};assembled={}
def txt(name,values,header):
    v=np.asarray(values);np.savetxt(O/name,v,fmt='%.17e',header=header)
    assert np.array_equal(np.loadtxt(O/name),v),name
    verified[name]={'shape':list(v.shape),'roundtrip_exact':True,'SHA256':m.track.sha(O/name)}

for case in ('hard','enge'):
    for obs in ('11m','1km'):
        key=case+'_'+obs;d=np.load(O/f'assembled_{key}.npz');assembled[key]=d
        lam=d['lambda_nm'];om=d['omega'];W=d['energy_xy'];old=d['matlab_original_total'];ar=d['matlab_rad_xyz'];at=d['matlab_totalfield_xyz']
        xy=ar[:2].sum(0);total=ar.sum(0);fullxy=at[:2].sum(0);an=d['analytic_total'];wc=d['energy_coarser']
        values=np.column_stack([lam,om,W,d['energy_components'][0],d['energy_components'][1],old,xy,total,ar[2],fullxy,
            d['matlab_velocity_xyz'][:2].sum(0),d['matlab_interference_xyz'][:2].sum(0),an,W-old,W-xy,W-fullxy,
            W/an-1,W/xy-1,W/old-1,W/fullxy-1,wc,W-wc,W/(m.HBAR*om)])
        obsxyz=cfg['observers'][obs]['position_lab_m'];R=cfg['observers'][obs]['R_source_m']
        txt(f'SRW_tangent_{key}.txt',values,
            f'240MeV single electron; {case}; observer_lab_m={obsxyz}; R_source={R:.17g}m. Native SRW xy, lab frame, compensated gamma/B.\n'
            'No fitted shift, scaling, smoothing, or Ez reconstruction. Energies J/(rad/s)/sr, photon photons/(rad/s)/sr.\n'
            'Columns: lambda_nm omega_rad_s SRW_xy SRW_x SRW_y original_MAGNETO_total MATLAB_fine_radiation_xy MATLAB_fine_radiation_total MATLAB_fine_radiation_z MATLAB_fine_radiation_plus_velocity_xy MATLAB_fine_velocity_xy MATLAB_fine_interference_xy analytic_total SRW_minus_original_total SRW_minus_fine_xy SRW_minus_LW_xy SRW_div_analytic_minus1 SRW_div_fine_xy_minus1 SRW_div_original_total_minus1 SRW_div_LW_xy_minus1 SRW_ds3p75um_xy SRW_ds1p875_minus_ds3p75 SRW_photon_xy')
        ld=d['dense_lambda_nm'];od=2*np.pi*m.C/(ld*1e-9);wd=d['dense_energy_xy'];rd=d['dense_matlab_rad_xyz'];td=d['dense_matlab_totalfield_xyz']
        txt(f'dense_{key}.txt',np.column_stack([ld,od,wd,rd[:2].sum(0),rd.sum(0),td[:2].sum(0),d['dense_analytic_total']]),
            'All9999 nodes directly calculated by SRW and new MATLAB reference.\nColumns: lambda_nm omega_rad_s SRW_xy MATLAB_radiation_xy MATLAB_radiation_total MATLAB_radiation_plus_velocity_xy analytic_total')
        txt(f'photon_{key}.txt',np.column_stack([lam,W/(m.HBAR*om),old/(m.HBAR*om),xy/(m.HBAR*om),fullxy/(m.HBAR*om),an/(m.HBAR*om)]),
            'Convert using hbar*omega, not input photon eV.\nColumns: lambda_nm SRW_xy original_MAGNETO_total MATLAB_fine_radiation_xy MATLAB_fine_radiation_plus_velocity_xy analytic_total')
        txt(f'integration_{key}.txt',np.column_stack([lam,W,wc,W-wc,W/wc-1]),
            'Full5000-node integral check, same trajectory and magnetic field.\nColumns: lambda_nm SRW_ds1p875um SRW_ds3p75um fine_minus_coarse fine_div_coarse_minus1')
        reconstruct=np.interp(np.log(ld[1::2]),np.log(lam),W)
        txt(f'wavelength_midpoints_{key}.txt',np.column_stack([ld[1::2],wd[1::2],reconstruct,wd[1::2]-reconstruct]),
            'Interpolation ONLY for this grid-resolution diagnostic; never used in formal spectra.\nColumns: lambda_mid_nm direct_SRW_xy loglambda_linear_reconstruction direct_minus_reconstruction')
        grid='refinement_'+key;ids=np.load(O/f'indices_{grid}.npy');lr=lam[ids]
        _,_,pf,_=r.collect(case,obs,grid,.9375e-6,3.75);_,_,pt,_=r.collect(case,obs,grid,1.875e-6,1.875)
        wf=pf.sum(0);wt=pt.sum(0)
        refinement[key]={'nodes':len(ids),'integration1p875_vs0p9375':r.metrics(W[ids],wf,lr),
            'trajectory_input3p75_vs1p875':r.metrics(W[ids],wt,lr)}
        txt(f'selected_refinement_{key}.txt',np.column_stack([lr,W[ids],wf,wt,W[ids]/wf-1,W[ids]/wt-1]),
            'Points selected after full-grid integration/benchmark comparisons plus evenly spread/end nodes. Hard trajectories retain the same fine force solution; Enge compares independently tracked original sources.\nColumns: lambda_nm main_ds1p875 finer_ds0p9375 finer_trajectory_input1p875 main_div_finer_integral_minus1 main_div_finer_trajectory_minus1')
        # These inputs were independently saved and exactly round-trip checked by their writers.
        for filename in (f'pilot_{key}.txt',f'band_grid_{key}.txt',f'velocity_incoming_window_bound_{key}.txt'):
            a=np.loadtxt(O/filename);assert np.all(np.isfinite(a))
            verified[filename]={'shape':list(a.shape),'roundtrip_exact':True,'SHA256':m.track.sha(O/filename),'checked_by':'pilot_review.py/results.py/velocity_window_bound.py'}

(O/'refinement_summary.json').write_text(json.dumps(refinement,indent=2))

# Export the exact production hard trajectory. Full fine force outputs remain in NPZ.
hard=np.load(O/'hard_b100000_ct3.75_edge0.125.npz')
cols=['ct_m',*m.track.COMPONENTS]
txt('hard_trajectory_radiation_input_3p75um.txt',np.column_stack([hard[k] for k in cols]),
    'Hard b1=1e5 model; force integration0.9375um, exact retained radiation input nodes3.75um. Native SRW track, gamma/B compensated.\nColumns: '+' '.join(cols))
field=np.load(O/'field_b100000_edge0.125.npz');rows=[]
for j in range(3):
    z=field[f'z{j}'];zz,yy,xx=np.meshgrid(z+1,field['y'],field['x'],indexing='ij')
    rows.append(np.column_stack([np.full(xx.size,j),xx.ravel(),yy.ravel(),zz.ravel(),field[f'B{j}'].reshape(-1,3)]))
txt('hard_fieldmap_3D_compensated.txt',np.concatenate(rows),
    'Three adjacent uniform SRW submaps, full GPT formula at all samples, shared boundary nodes retained. Field already compensated: DO NOT rescale again.\nColumns: submap_id x_WCS_m y_WCS_m z_WCS_m Bx_T By_T Bz_T')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
plt.rcParams.update({'font.size':10,'axes.formatter.use_mathtext':True})
def style(ax):
    ax.set_xscale('log');ax.set_xlim(15,150);ax.set_xticks([15,20,30,40,50,70,100,150]);ax.set_xticklabels(['15','20','30','40','50','70','100','150'])
    ax.xaxis.set_minor_formatter(NullFormatter());ax.grid(alpha=.25)
def save(fig,name):
    fig.savefig(O/(name+'.png'),dpi=190);fig.savefig(O/(name+'.svg'));plt.close(fig)
colors={'hard':'#d95f02','enge':'#7b3294'}
for obs in ('11m','1km'):
    fig,ax=plt.subplots(3,1,figsize=(12,11),sharex=True,layout='constrained')
    an=assembled['hard_'+obs]['analytic_total'];lam=assembled['hard_'+obs]['lambda_nm']
    ax[0].loglog(lam,an,'k',lw=1.6,label='Analytic: infinite uniform arc, total')
    for case in ('hard','enge'):
        d=assembled[case+'_'+obs];old=d['matlab_original_total'];W=d['energy_xy'];xy=d['matlab_rad_xyz'][:2].sum(0)
        ax[0].loglog(lam,old,ls='--' if case=='hard' else ':',lw=1.3,color='#1f77b4' if case=='hard' else '#238b45',label=f'Original MAGNETO {case}: total')
        ax[0].loglog(lam,W,color=colors[case],lw=.85,label=f'Native SRW {case}: x+y')
        ax[1].plot(lam,(W/an-1)*100,color=colors[case],lw=.8,label=f'SRW {case} xy / analytic total - 1')
        ax[1].plot(lam,(old/an-1)*100,color='#1f77b4' if case=='hard' else '#238b45',lw=.6,label=f'MAGNETO {case} total / analytic - 1')
        ax[2].plot(lam,(W/old-1)*100,color=colors[case],lw=.8,label=f'SRW {case} xy / original MAGNETO total - 1')
        ax[2].plot(lam,(W/xy-1)*100,color=colors[case],ls='--',lw=.8,label=f'SRW {case} xy / fine MATLAB same xy - 1')
    ax[0].set_ylabel('Energy density [J/(rad/s)/sr]');ax[0].legend(fontsize=8,loc='lower right')
    ax[1].set_ylabel('Difference from analytic [%]');ax[1].legend(fontsize=8)
    ax[2].set_ylabel('Difference from MATLAB [%]');ax[2].set_xlabel('Wavelength [nm]');ax[2].legend(fontsize=8)
    for aa in ax:style(aa)
    theta=cfg['tangent_theta_rad']*1e3;R=cfg['observers'][obs]['R_source_m']
    fig.suptitle(f'240 MeV dipole, midpoint tangent {theta:.6f} mrad, observer {obs}\nB0=0.65 T, L=0.1 m; per steradian from fixed tangent point, R_source={R:.6f} m')
    save(fig,f'compare_tangent_{obs}')

    fig,ax=plt.subplots(2,1,figsize=(12,7),sharex=True,layout='constrained')
    for case in ('hard','enge'):
        key=case+'_'+obs;d=assembled[key];W=d['energy_xy'];wc=d['energy_coarser']
        ids=np.load(O/f'indices_refinement_{key}.npy');_,_,pf,_=r.collect(case,obs,'refinement_'+key,.9375e-6,3.75)
        ax[0].plot(lam,(W/wc-1)*1e6,lw=.7,color=colors[case],label=f'{case}: full5000 nodes, ds1.875 / ds3.75 - 1')
        ax[0].scatter(lam[ids],(W[ids]/pf.sum(0)-1)*1e6,s=10,color=colors[case],label=f'{case}: selected, ds1.875 / ds0.9375 - 1')
        ld=d['dense_lambda_nm'];wd=d['dense_energy_xy'];interp=np.interp(np.log(ld[1::2]),np.log(lam),W)
        ax[1].plot(ld[1::2],(wd[1::2]-interp)/max(wd)*1e6,lw=.65,color=colors[case],label=f'{case}: direct midpoint - reconstructed / peak')
    ax[0].set_ylabel('Integral difference [ppm]');ax[0].legend(fontsize=8)
    ax[1].set_ylabel('Grid reconstruction error / peak [ppm]');ax[1].set_xlabel('Wavelength [nm]');ax[1].legend(fontsize=8)
    for aa in ax:style(aa)
    fig.suptitle(f'240 MeV dipole tangent, {obs}: numerical checks\nInterpolation is only a grid diagnostic; formal spectra are direct native SRW')
    save(fig,f'convergence_tangent_{obs}')

fig,ax=plt.subplots(2,1,figsize=(11,7),layout='constrained')
g=np.load(O/'GPT_hard.npz');fine=np.load(O/'hard_b100000_ct0.9375_edge0.125.npz')
inside=fine['ct_m']>=g['ct_m'][0]
ax[0].plot(g['z_m'],g['y_m']*1e3,label='GPT hard reference');ax[0].plot(fine['z_m'][inside],fine['y_m'][inside]*1e3,'--',label='Independent SRW, force ct0.9375um')
ax[0].set_ylabel('Vertical orbit [mm]');ax[0].legend(fontsize=8)
ax[1].plot(g['z_m'],(np.interp(g['ct_m'],fine['ct_m'],fine['y_m'])-g['y_m'])*1e12,lw=.8)
ax[1].set_ylabel('SRW - GPT y [pm]')
for aa in ax:aa.set_xlabel('Global z [m]');aa.grid(alpha=.25)
fig.suptitle('New hard dipole independent SRW trajectory, b1=1e5 /m\nNo fitted field, position, angle or time shift')
save(fig,'hard_dipole_tracking')

(O/'txt_verified.json').write_text(json.dumps(verified,indent=2))
provenance={'protected_inputs_unchanged':all(m.track.sha(p)==h for p,h in cfg['protected_SHA256'].items()),
    'script_SHA256':{str(p):m.track.sha(p) for p in O.glob('*.py')},
    'MATLAB_script_SHA256':{str(p):m.track.sha(p) for p in O.glob('*.m')},
    'trajectory_and_field_SHA256':{str(p):m.track.sha(p) for p in [*O.glob('hard_*.npz'),*O.glob('field_*.npz'),*O.glob('force_refined_*.npz')]},
    'native_binary_SHA256':m.track.sha(m.track.SRW/'env/python/srwpy/srwlpy.pyd')}
assert provenance['protected_inputs_unchanged'];(O/'provenance.json').write_text(json.dumps(provenance,indent=2))
print(json.dumps({'TXT_files':len(verified),'selected_refinement':refinement},indent=2),flush=True)
