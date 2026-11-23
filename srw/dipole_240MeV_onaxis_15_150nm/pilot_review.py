from pathlib import Path
import sys,json
sys.dont_write_bytecode=True
import calculate as m
import numpy as np
from scipy.io import loadmat
O=m.O
def spectrum(ct=1.875,fs=2.5,ds=.46875e-6,term=1):
    folder=O/'blocks'/f'pilot_ct{ct:g}_field{fs:g}_ds{ds:g}_full_t{term}'
    if (folder/'000000_000090.npz').exists():
        a=np.load(folder/'000000_000090.npz');return a['energy_components'].sum(axis=0)
    return np.concatenate([np.load(folder/name)['energy_components'].sum(axis=0) for name in ('000000_000044.npz','000044_000090.npz')])
lam=np.load(O/'lambda_pilot.npy');ids=np.load(O/'pilot_indices.npy')
ref=np.load(O/'reference_original.npz');old=ref['energy_xy'][ids]
fine=spectrum(ds=.1171875e-6);mat=loadmat(O/'matlab_pilot_ct1.875um.mat',simplify_cells=True)['result'];new=mat['W_total_field'][:2].sum(axis=0)
def check(a,b):return {'max_pointwise_relative':float(max(abs(a/b-1))), 'max_difference_over_peak':float(max(abs(a-b))/max(b)),
    'at15nm_relative':float(a[0]/b[0]-1),'at150nm_relative':float(a[-1]/b[-1]-1)}
checks={'integration_vs_ds0p1171875':{f'{ds:g}':check(spectrum(ds=ds),fine) for ds in (3.75e-6,1.875e-6,.9375e-6,.46875e-6,.234375e-6)},
    'trajectory_ct3p75_vs1p875':check(spectrum(ct=3.75),spectrum()),'trajectory_ct1p875_vs0p9375':check(spectrum(),spectrum(ct=.9375)),
    'field2p5_vs1p25':check(spectrum(ct=.9375),spectrum(ct=.9375,fs=1.25)),
    'native_vs_original_matlab':check(fine,old),'native_vs_velocity_inclusive_reference':check(fine,new),
    'termination0_vs1':check(spectrum(term=0),spectrum()),
    'velocity_inclusive_vs_original':check(new,old),
    'native_implementation':'The GPT-output-only sub-window with original ct=0 initial metadata is rejected by native SRW because initial conditions lie outside the trajectory domain. No spectrum was produced for that diagnostic. All full original ct production/pilot calls succeed.',
    'production_ct_um':1.875,'production_field_um':2.5,'production_ds_m':.234375e-6}
(O/'pilot_review.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
data=np.column_stack([lam,old,new,*[spectrum(ds=d) for d in (3.75e-6,1.875e-6,.9375e-6,.46875e-6,.234375e-6,.1171875e-6)],spectrum(term=0)])
np.savetxt(O/'pilot_comparison.txt',data,fmt='%.17e',
    header='Columns: lambda_nm original_MATLAB_xy new_LW_radiation_plus_velocity_xy SRW_ds3p75um SRW_ds1p875um SRW_ds0p9375um SRW_ds0p46875um SRW_ds0p234375um SRW_ds0p1171875um SRW_ds0p46875um_term0')
assert np.array_equal(np.loadtxt(O/'pilot_comparison.txt'),data)
