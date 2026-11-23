"""Exact upper bound on omitted incoming uniform-motion velocity-field Fourier term.
No asymptotic approximation is used: integration by parts plus monotone |Ey(t_obs)|.
This bounds a reference-window difference; it does not modify any formal spectrum.
"""
import sys,json
sys.dont_write_bytecode=True
import model as m
import numpy as np
from scipy.io import loadmat
cfg=json.loads((m.O/'manifest.json').read_text());eps=8.85418781e-12;q=-1.6021892e-19
gamma=cfg['gamma_GPT'];beta=np.sqrt(1-gamma**-2);out={}
for case in ('hard','enge'):
    for obsname in ('11m','1km'):
        obs=np.asarray(cfg['observers'][obsname]['position_lab_m']);R=np.linalg.norm(obs);n=obs/R;D=1-n[2]*beta
        Ey0=q/(4*np.pi*eps)*n[1]/(gamma**2*R**2*D**3)
        # For the exact constant-velocity field, transverse |Ey(tobs)| decreases
        # once the virtual instantaneous source is downstream of the observation z.
        # These are algebraic conditions for the LW field of the known incoming segment,
        # not an assumption about the actual particle's future path after the bend.
        assert beta*R>obs[2]
        a=loadmat(m.O/f'MATLAB_baseline_{case}_{obsname}_ct1.875um.mat',simplify_cells=True)['result']
        omega=np.asarray(a['omega']);F=a['Ew_radiation_lab'][1]+a['Ew_velocity_lab'][1]
        bound=2*abs(Ey0)/omega
        ratio=bound/abs(F);impact=2*ratio+ratio**2
        name=f'velocity_incoming_window_bound_{case}_{obsname}.txt'
        data=np.column_stack([a['lambda_nm'],bound,abs(F),impact])
        np.savetxt(m.O/name,data,fmt='%.17e',header='Rigorous incoming-uniform-segment bound only; no values added to formal spectrum.\nColumns: lambda_nm max_abs_missing_Fy_Vs_per_m abs_existing_Fy_Vs_per_m max_relative_xy_energy_change_bound')
        assert np.array_equal(np.loadtxt(m.O/name),data)
        out[case+'_'+obsname]={'abs_Ey_at_first_arrival_V_per_m':float(abs(Ey0)),
            'monotonicity_condition_beta_normobs_minus_obs_z_m':float(beta*R-obs[2]),
            'max_relative_xy_spectrum_change_bound':float(max(impact)),
            'method':'For known original field-free incoming straight segment, |Ey(tobs)| is monotone decreasing. Integration by parts gives |integral Ey exp(-i omega tobs) dtobs| <= 2|Ey(first arrival)|/omega. Additive field perturbation changes |F|^2 by at most 2a/|F|+(a/|F|)^2.',
            'no_new_physics':'Uses the original GPT initial condition and its existing field-free ct0..first-record interval, no external drift, no physical future straight-line assumption, no approximation or spectrum correction.'}
(m.O/'velocity_window_bound.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
