"""Summarize stage-one results without changing the native SRW calculation."""
from pathlib import Path
import json
import sys
import os
import hashlib
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT.parent/'.srw_dipole_deps'))
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'.matplotlib'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

inputs=json.loads((OUT/'inputs.json').read_text())
extra=[OUT/'track.py', OUT/'finalize_tracking.py',
    OUT.parent/'SRW-master/SRW-master/cpp/src/core/srebmdat.h',
    OUT.parent/'SRW-master/SRW-master/cpp/src/lib/srwlib.cpp']
source=Path(__file__).resolve().parents[1]/'gpt_sims'/'single_electron'
for folder,inp in [('E7GeV_und15p_K1','single_electron_und15p.in'),('E7GeV_dipole_z1m','single_electron_dipole.in')]:
    extra += [source/folder/'stats.gdf',source/folder/'result.gdf',source/folder/inp]
for path in extra:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    inputs['inputs'][str(path)]=h.hexdigest()
(OUT/'inputs.json').write_text(json.dumps(inputs,indent=2))
summary=json.loads((OUT/'tracking_summary.json').read_text())
gamma=inputs['reference']['dipole']['gamma']
rest_interface=.51099890221
rest_internal=.511003414
internal_gamma=gamma*rest_interface/rest_internal
gpt=inputs['reference']['dipole']
ratio=(1.602176462e-19/9.10938188e-31/internal_gamma)/(abs(gpt['q'])/gpt['m']/gamma)
diagnostics={'input_gamma':gamma,'SRW_internal_gamma_from_source':internal_gamma,
    'SRW_interface_rest_energy_MeV':rest_interface,'SRW_internal_rest_energy_MeV':rest_internal,
    'relative_gamma_shift':internal_gamma/gamma-1,
    'SRW_over_GPT_magnetic_acceleration_coefficient':ratio,
    'explanation':'Source-based prediction only; no gamma, mass, charge, field or output correction applied.'}
fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
for col,case in enumerate(('undulator','dipole')):
    suffix='' if case=='undulator' else '_field15um'
    a=dict(np.load(OUT/f'{case}_ct15um{suffix}.npz'))
    ref=dict(np.load(OUT/f'{case}_gpt_reference.npz'))
    coord='x' if case=='undulator' else 'y'
    thin=4
    z=ref['avgz']
    dp=a[coord+'_m'][::thin]-ref['avg'+coord]
    dv=a['beta_'+coord][::thin]-ref['avgB'+coord]
    axes[0,col].plot(z,dp*1e9)
    axes[1,col].plot(z,dv*1e6)
    axes[0,col].set_title(case.capitalize())
    axes[0,col].set_ylabel(f'SRW - GPT: {coord} [nm]')
    axes[1,col].set_ylabel(r'SRW - GPT: $\beta_'+coord+r'$ [$\mu$rad]')
    for ax in axes[:,col]:
        ax.grid(alpha=.25); ax.set_xlabel('Global z [m]')
    final_angle=float(np.arctan2(a['beta_'+coord][-1],a['beta_z'][-1]))
    reference_angle=float(np.arctan2(ref['avgB'+coord][-1],ref['avgBz'][-1]))
    diagnostics[case]={'final_angle_rad':final_angle,'GPT_final_angle_rad':reference_angle,
        'final_angle_difference_rad':final_angle-reference_angle,
        'max_abs_transverse_position_m':float(np.max(np.abs(a[coord+'_m']))),
        'max_abs_transverse_beta':float(np.max(np.abs(a['beta_'+coord]))),
        'max_abs_transverse_beta_difference_over_peak':float(np.max(np.abs(dv))/np.max(np.abs(ref['avgB'+coord])))}
    if case=='dipole':
        r=np.array([a['x_m'][-1],a['y_m'][-1],a['z_m'][-1]])
        u=np.array([a['beta_x'][-1],a['beta_y'][-1],a['beta_z'][-1]])
        u/=np.linalg.norm(u)
        s=-r.dot(u)+np.sqrt(r.dot(u)**2+11**2-r.dot(r))
        hit=r+s*u
        diagnostics[case]['straight_exit_ray_sphere_hit_m']=hit.tolist()
fig.suptitle('Native SRW versus GPT at identical original output times (no fitted correction)')
fig.savefig(OUT/'srw_minus_gpt.png',dpi=170);plt.close(fig)
(OUT/'constants_diagnostic.json').write_text(json.dumps(diagnostics,indent=2))

lines=['# AAC benchmark — stage 1: independent SRW particle tracking','',
'Only particle tracking has been performed. Radiation calculations have not started.', '',
'## Physical inputs and numerical representation','',
'- One electron starts at (0,0,0), parallel to +z, in each independent case.',
f'- GPT kinetic energy: 7 GeV. Exact gamma read from GPT: {gamma:.12f}.',
'- Undulator: complete original 31 × 21 × 3001 magnetic map; mm converted to m; every field value preserved. Active +90° rotation about z: (x,y,z)=(-y_local,x_local,z_local+1); B=(-By_local,Bx_local,Bz_local). Electron wiggles horizontally.',
'- Magnetic interpolation is SRW trilinear (interp=1), matching GPT. Original grid coordinates are supplied explicitly; no sinusoidal fit, field normalization, fringe removal, smoothing or field-grid refinement is applied to the undulator.',
'- Dipole: the user-approved table is evaluated from the installed GPT rectmagnet_sim source, not SRW default dipole fringes. Original B0=0.059312008063175437 T, b1=40/m, b2=dl=0, half-length=0.5 m, hard truncation 0.25 m beyond either edge, and width 0.2 m are retained.',
'- Dipole table coordinates: global x = (-1,0,+1) mm, y = (-0.1,0,+0.1) m. Formula includes transverse Bz and x dependence; it is independent of y inside the original width. The exact reference orbit remains on the x=0 node, so no transverse reduction or averaging is used on that orbit. This table is for this central single-electron orbit, not a validation for off-axis particle ensembles.',
'- SRW CalcPartTraj with manual fixed-step Runge–Kutta tracking ([1]). ct=c*t is the independent variable. Baseline samples match original GPT times; subsequent grids are nested at half and quarter of that step. Tiny accumulated GPT output-time rounding is at sub-picometre ct scale.',
'- Step tests: ct=60,30,15 um. Dipole field tables: dz=60,30,15 um, tested at fixed ct=15 um. Tracking-step tests keep dz=15 um fixed.',
'- Native SRW charge, electron mass, and relE0=1 are retained. No fitting or compensating adjustment has been made. No bunch model, collective effects, field-free extrapolation, or custom particle integrator is introduced.', '',
'## Important native SRW constants discrepancy','',
f'The C++ interface uses electron rest energy {rest_interface} MeV to convert the supplied gamma to energy; srTEbmDat::SetupGamma uses {rest_internal} MeV to convert it back. Therefore the effective internal gamma is {internal_gamma:.12f}, a relative shift of {(internal_gamma/gamma-1)*1e6:.6f} ppm. The saved initial particle metadata still contains the supplied gamma; it must not be mistaken for the internal integration gamma.', '',
f'Together with the differing native SRW/GPT q/m constants, the magnetic acceleration coefficient q/(gamma*m*c) is predicted to be larger in SRW by {(ratio-1)*1e6:.6f} ppm. The observed dipole bend discrepancy has the same size. This is a source-based explanation, not a correction applied to these outputs.', '',
'The actual native speed-norm residual relative to the supplied gamma is about 9.43e-14 in beta squared, consistent with the internal gamma difference. Constant beta norm does not independently validate the trajectory or subsequent radiation phase accuracy.', '',
'## Results and measured changes','',
'All differences below are maximum absolute differences along the common trajectory samples. They are not relative errors at trajectory zero crossings. No pass/fail convergence tolerance has been imposed.', '',
'| Case | Max transverse position difference from GPT | Max transverse beta difference from GPT | Final angle SRW / GPT |',
'|---|---:|---:|---:|']
for case in ('undulator','dipole'):
    coord='x' if case=='undulator' else 'y'
    d=summary[case]['finest_vs_GPT_at_original_times'];v=diagnostics[case]
    lines.append(f"| {case} | {d[coord+'_m']*1e9:.8g} nm | {d['beta_'+coord]:.8g} | {v['final_angle_rad']*1e3:.10g} / {v['GPT_final_angle_rad']*1e3:.10g} mrad |")
lines+=['','| Case | ct comparison | Max transverse orbit change | Max transverse beta change |',
'|---|---|---:|---:|']
for case in ('undulator','dipole'):
    coord='x' if case=='undulator' else 'y'
    for row in summary[case]['tracking_convergence_fixed_field'][:-1]:
        d=row['vs_finest_max_absolute']
        lines.append(f"| {case} | {row['ct_step_um']:.0f} vs 15 um | {d[coord+'_m']*1e9:.8g} nm | {d['beta_'+coord]:.8g} |")
lines+=['','| Dipole field table comparison, fixed ct=15 um | Max y change | Max beta_y change |',
'|---|---:|---:|']
for row in summary['dipole']['field_convergence_fixed_ct15um'][:-1]:
    d=row['vs_field15um_max_absolute']
    lines.append(f"| dz={row['field_step_um']:g} vs 15 um | {d['y_m']*1e9:.8g} nm | {d['beta_y']:.8g} |")
lines+=['',
'The undulator position ranges from about -0.4022 to +0.4948 um; y and beta_y remain exactly zero. Dipole x and beta_x remain exactly zero; final y is about 2.0318765 mm. The straight exit ray reaches the R=11 m sphere at y='+f"{diagnostics['dipole']['straight_exit_ray_sphere_hit_m'][1]*1e3:.8f} mm.",
'',
'Refinement differences are small but not strictly monotonic, particularly for the undulator with piecewise-linear field interpolation and truncated map boundaries. These comparisons demonstrate the measured stability at these three grids; they are not an extrapolated rigorous bound. Radiation convergence still needs separate frequency-, angle-, trajectory-, and integration-precision checks.', '',
'![Undulator trajectory](undulator_trajectory.png)','',
'![Dipole trajectory](dipole_trajectory.png)','',
'![Native SRW minus GPT](srw_minus_gpt.png)','',
'## Reproduction and files','',
'Use Python 3.12. SRW is the existing local srwlpy.pyd; binary and input SHA256 hashes are recorded in inputs.json. The original D: files were read, not modified.', '',
'Run track.py prepare; track.py run --case undulator --factor 1 (then 2 and 4); track.py run --case dipole --factor 1 --fieldstep 15 (then factors 2 and 4); at factor 4 additionally run --fieldstep 60 and 30; then track.py report; finalize_tracking.py. Runs can be wrapped by ../srw_dipole_simulation/guarded_run.py for process limits.', '',
'- *_trajectory.csv: finest computed independent SRW orbit, velocities, fields, and ct; full double-precision numbers.',
'- *_ct*.npz and *.json: every grid and its timing/settings; *_gpt_reference.npz: GPT reference arrays.',
'- tracking_summary.json: absolute comparison metrics and all numerical checks.',
'- constants_diagnostic.json: native constants mapping, bend angles and mirror intersection.',
'- track.py and finalize_tracking.py: complete reproducible calculation and report generation.', '',
'## Pending before radiation','',
'Review these tracking results. Decide explicitly whether subsequent comparisons should preserve native SRW constants or compensate the documented interface gamma conversion to match physical gamma. Such compensation has not been applied. No radiation computation is authorized to continue automatically past this review checkpoint.']
(OUT/'TRACKING_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(diagnostics,indent=2))
