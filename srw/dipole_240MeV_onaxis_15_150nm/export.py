"""Export exact unsmoothed results and publication figures with numerical limitations visible."""
import sys,json
from pathlib import Path
sys.dont_write_bytecode=True
import calculate as m
import analysis as a
import numpy as np
O=m.O
s=json.loads((O/'summary.json').read_text());x=np.load(O/'assembled.npz')
lam=x['lambda_nm'];om=x['omega'];w=x['energy_xy'];old=x['matlab_original_xy'];new=x['matlab_total_field_xy'];coarse=x['energy_coarser']
ri=np.load(O/'refinement_indices.npy');lr=lam[ri]
_,f,sec=a.collect('refinement',.1171875e-6);wr=(abs(f)**2*m.HBAR*1e9*121).sum(0)
_,ft,_=a.collect('refinement',.234375e-6,.9375,2.5);wt=(abs(ft)**2*m.HBAR*1e9*121).sum(0)
_,ff,_=a.collect('refinement',.234375e-6,.9375,1.25);wf=(abs(ff)**2*m.HBAR*1e9*121).sum(0)
refined={'nodes':len(ri),'integration_ds0p234375_vs0p1171875':a.metrics(w[ri],wr,lr),
    'trajectory_ct1p875_vs0p9375_same_ds0p234375':a.metrics(w[ri],wt,lr),
    'field2p5_vs1p25_same_ct0p9375_same_ds0p234375':a.metrics(wt,wf,lr)}
(O/'refinement_summary.json').write_text(json.dumps(refined,indent=2));print(json.dumps(refined,indent=2),flush=True)
verified={}
def txt(name,values,header):
    values=np.asarray(values);np.savetxt(O/name,values,fmt='%.17e',header=header)
    check=np.loadtxt(O/name);assert np.array_equal(check,values),name
    verified[name]={'shape':list(values.shape),'SHA256':m.track.sha(O/name),'roundtrip_exact':True}
pol=x['energy_components'];totalxyz=x['matlab_total_field_xyz']
values=np.column_stack([lam,om,w,pol[0],pol[1],old,x['matlab_original_total'],x['matlab_original_z'],
    x['matlab_radiation_xyz'][:2].sum(0),x['matlab_velocity_xyz'][:2].sum(0),new,totalxyz.sum(0),x['matlab_interference_xyz'][:2].sum(0),
    w-old,w-new,w/(m.HBAR*om),coarse,w-coarse,w/old-1,w/new-1])
txt('dipole240_onaxis_spectrum.txt',values,
    '240MeV single electron, on axis (0,0,11)m. Direct native SRW; no fitted shift, scaling, smoothing or point replacement. Energy per angular frequency, NOT per wavelength.\n'
    'Energy columns J/(rad/s)/sr, photon column photons/(rad/s)/sr.\n'
    'Columns: lambda_nm omega_rad_s SRW_xy SRW_x SRW_y MATLAB_original_xy MATLAB_original_total MATLAB_original_z LW_radiation_xy LW_velocity_xy LW_radiation_plus_velocity_xy LW_radiation_plus_velocity_xyz LW_interference_xy SRW_minus_MATLAB_xy SRW_minus_LW_xy SRW_photon_xy SRW_ds0p46875um_xy SRW_ds0p234375_minus_ds0p46875 SRW_div_MATLAB_xy_minus1 SRW_div_LW_xy_minus1')
ld=x['dense_lambda_nm'];wd=x['dense_energy_xy'];nd=x['dense_matlab_total_field_xy'];od=2*np.pi*m.C/(ld*1e-9)
txt('dipole240_dense_grid.txt',np.column_stack([ld,od,wd,nd,wd-nd,wd/nd-1]),
    'All9999 wavelengths independently computed by both native SRW and added MATLAB exact LW reference.\nColumns: lambda_nm omega_rad_s SRW_xy LW_radiation_plus_velocity_xy SRW_minus_LW_xy SRW_div_LW_xy_minus1')
txt('dipole240_photon_spectrum.txt',np.column_stack([lam,w/(m.HBAR*om),old/(m.HBAR*om),new/(m.HBAR*om)]),
    'Photon conversion uses hbar*omega, no input-eV conversion.\nColumns: lambda_nm SRW_photon_xy MATLAB_original_photon_xy LW_radiation_plus_velocity_photon_xy')
txt('difference_comparisons.txt',np.column_stack([lam,w-old,w-x['matlab_original_total'],w-new,
    w/old-1,w/x['matlab_original_total']-1,w/new-1,x['matlab_original_z'],totalxyz[2]]),
    'Native SRW xy retained; original total includes z. No longitudinal SRW reconstruction.\nColumns: lambda_nm SRW_minus_original_MATLAB_xy SRW_minus_original_MATLAB_total SRW_minus_added_LW_xy SRW_div_original_xy_minus1 SRW_div_original_total_minus1 SRW_div_added_LW_xy_minus1 original_MATLAB_z added_LW_z')
txt('integration_convergence.txt',np.column_stack([lam,w,coarse,w-coarse,w/coarse-1]),
    'Full original grid, same trajectory and field; no discarded weak nodes.\nColumns: lambda_nm SRW_ds0p234375um SRW_ds0p46875um fine_minus_coarse fine_div_coarse_minus1')
txt('selected_refinement.txt',np.column_stack([lr,w[ri],wr,wt,wf,w[ri]/wr-1,w[ri]/wt-1,wt/wf-1]),
    'Selected only after examining complete-grid differences: even-spaced nodes and largest integration/benchmark discrepancies.\nColumns: lambda_nm main_ds0p234375 finer_ds0p1171875 ct0p9375_ds0p234375 field1p25_ct0p9375_ds0p234375 main_div_finer_minus1 main_div_finer_ct_minus1 field2p5_div_field1p25_minus1')
txt('wavelength_midpoint_check.txt',np.column_stack([ld[1::2],wd[1::2],np.interp(np.log(ld[1::2]),np.log(lam),w),
    wd[1::2]-np.interp(np.log(ld[1::2]),np.log(lam),w)]),
    'Interpolation ONLY for this grid-resolution diagnostic, not formal output.\nColumns: lambda_mid_nm directly_computed_SRW_xy loglambda_linear_reconstruction_of_original_grid direct_minus_reconstruction')
pilot=json.loads((O/'pilot_review.json').read_text())
checks=np.loadtxt(O/'pilot_comparison.txt');verified['pilot_comparison.txt']={'shape':list(checks.shape),'SHA256':m.track.sha(O/'pilot_comparison.txt'),'roundtrip_exact':True,'verified_by':'pilot_review.py'}
for filename in ('band_grid_checks.txt',):
    v=np.loadtxt(O/filename);assert np.all(np.isfinite(v));verified[filename]={'shape':list(v.shape),'SHA256':m.track.sha(O/filename),'roundtrip_exact':True,'verified_by':'analysis.py'}
(O/'txt_verified.json').write_text(json.dumps(verified,indent=2))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.formatter.use_mathtext':True})
colors={'srw':'#d95f02','old':'#1f77b4','lw':'#202020'}
def finish(fig,name):
    fig.savefig(O/(name+'.png'),dpi=190);fig.savefig(O/(name+'.svg'));plt.close(fig)
def axes_style(ax):
    ax.set_xscale('log');ax.set_xlim(15,150);ax.set_xticks([15,20,30,50,70,100,150]);ax.set_xticklabels(['15','20','30','50','70','100','150']);ax.grid(alpha=.25)

fig,ax=plt.subplots(3,1,figsize=(11,11),sharex=True,layout='constrained')
ax[0].loglog(lam,old,color=colors['old'],lw=1.7,label='Original MATLAB: LW acceleration, x+y')
ax[0].loglog(lam,new,color=colors['lw'],lw=1.25,ls='--',label='Added MATLAB: LW acceleration + velocity, x+y')
ax[0].loglog(lam,w,color=colors['srw'],lw=.9,label='Native SRW: x+y')
ax[0].set_ylabel(r'$d^2W/(d\omega\,d\Omega)$ [J/(rad/s)/sr]');ax[0].legend(fontsize=8,loc='upper left')
ax[1].plot(lam,w-old,color=colors['old'],lw=.9,label='SRW - original MATLAB x+y')
ax[1].plot(lam,w-new,color=colors['lw'],lw=.9,label='SRW - acceleration+velocity reference x+y')
ax[1].set_ylabel('Signed spectral difference');ax[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax[1].legend(fontsize=8)
ax[2].plot(lam,(w/old-1)*100,color=colors['old'],lw=.65,label='SRW / original MATLAB - 1')
ax[2].plot(lam,(w/new-1)*100,color=colors['lw'],lw=.65,label='SRW / acceleration+velocity - 1')
ax[2].axhline(0,color='.5',lw=.6);ax[2].set_ylabel('Pointwise difference [%]');ax[2].set_xlabel('Wavelength [nm]');ax[2].legend(fontsize=8)
for aa in ax:axes_style(aa)
fig.suptitle('240 MeV short dipole: on-axis benchmark at (0, 0, 11) m\nB0=0.65 T, L=0.1 m, gap=10 mm; original 5000 wavelengths; no fitted correction')
finish(fig,'dipole240_SRW_MATLAB_comparison')

fig,ax=plt.subplots(2,1,figsize=(11,7),sharex=True,layout='constrained')
ax[0].plot(lam,(w/coarse-1)*100,lw=.7,label='Full5000 nodes: ds0.234375um / ds0.46875um - 1')
ax[0].scatter(lr,(w[ri]/wr-1)*100,s=12,color='#d95f02',label='Selected nodes: ds0.234375um / ds0.1171875um - 1',zorder=3)
ax[0].axhline(0,color='.5',lw=.6);ax[0].set_ylabel('Integration difference [%]');ax[0].legend(fontsize=8)
ax[1].plot(lr,(w[ri]/wt-1)*100,'.-',lw=.6,ms=3,label='ct1.875um / ct0.9375um - 1; same ds, field')
ax[1].plot(lr,(wt/wf-1)*100,'.-',lw=.6,ms=3,label='field2.5um / field1.25um - 1; same ds, ct')
ax[1].set_ylabel('Tracking / field difference [%]');ax[1].set_xlabel('Wavelength [nm]');ax[1].legend(fontsize=8)
for aa in ax:axes_style(aa)
fig.suptitle('Numerical checks: native radiation integral, trajectory sampling and magnetic field table\nShort-wave weak signals retain non-monotonic numerical sensitivity; no nodes removed')
finish(fig,'dipole240_convergence')

fig,ax=plt.subplots(2,1,figsize=(11,7),sharex=True,layout='constrained')
ax[0].plot(lam,old,color=colors['old'],label='Original MATLAB: acceleration only')
ax[0].plot(lam,new,'--',color=colors['lw'],label='Added MATLAB: acceleration + velocity')
ax[0].plot(lam,w,color=colors['srw'],lw=.8,label='Native SRW')
ax[0].set_ylabel('Energy density [J/(rad/s)/sr]');ax[0].legend(fontsize=8)
ax[1].plot(lam,(new/old-1)*100,color=colors['lw']);ax[1].set_ylabel('Velocity-field inclusion effect [%]');ax[1].set_xlabel('Wavelength [nm]')
for aa in ax:axes_style(aa)
fig.suptitle('240 MeV dipole: linear intensity scale and effect of including the LW velocity field')
finish(fig,'dipole240_linear_and_velocity_effect')

fig,ax=plt.subplots(2,1,figsize=(11,7),sharex=True,layout='constrained')
ax[0].plot(ld[1::2],(wd[1::2]-np.interp(np.log(ld[1::2]),np.log(lam),w))/max(wd),lw=.7)
ax[0].set_ylabel('Midpoint reconstruction error / peak')
ax[1].plot(lam,x['matlab_total_field_coarse']/new-1,lw=.7)
ax[1].set_ylabel('MATLAB ct3.75 / ct1.875 - 1');ax[1].set_xlabel('Wavelength [nm]')
for aa in ax:axes_style(aa)
fig.suptitle('Independent9999-node wavelength grid and added MATLAB trajectory-sampling check\nInterpolation is used only in the top diagnostic; all formal spectra are direct calculations')
finish(fig,'dipole240_grid_and_reference_checks')

provenance={'input_manifest_SHA256':m.track.sha(O/'manifest.json'),
    'all_radiation_protected_inputs_unchanged':all(m.track.sha(p)==h for p,h in json.loads((O/'manifest.json').read_text())['protected_SHA256'].items()),
    'all_tracking_original_inputs_unchanged':all(m.track.sha(p)==h for p,h in json.loads((O/'manifest.json').read_text())['tracking_original_SHA256'].items()),
    'scripts_SHA256':{str(p):m.track.sha(p) for p in O.glob('*.py')},
    'matlab_scripts_SHA256':{str(p):m.track.sha(p) for p in O.glob('*.m')},
    'native_binary_SHA256':m.track.sha(m.track.SRW/'env/python/srwpy/srwlpy.pyd')}
assert provenance['all_radiation_protected_inputs_unchanged'] and provenance['all_tracking_original_inputs_unchanged']
(O/'provenance.json').write_text(json.dumps(provenance,indent=2))
print('All spectra exported; PNG/SVG and exact TXT round trips complete.',flush=True)
