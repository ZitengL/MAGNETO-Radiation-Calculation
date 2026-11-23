import sys,json
sys.dont_write_bytecode=True
import model as m
import numpy as np
from scipy.io import loadmat
O=m.O;cfg=json.loads((O/'manifest.json').read_text());ref=np.load(O/'reference_original.npz')
ids=np.load(O/'pilot_indices.npy');lam=np.load(O/'lambda_pilot.npy')
def get(case,obs,ct,prec,edge=.125,method=0):
    label=f'pilot_{case}_{obs}_lab_ct{ct:g}_edge{edge:g}_b100000_m{method}_p{prec:g}'
    path=O/'blocks'/label/f'000000_{len(lam):06d}.npz'
    if not path.exists():return None
    return np.load(path)['energy_components'].sum(0)
def metrics(a,b):
    d=a/b-1;k=int(np.argmax(abs(d)))
    return {'max_abs_relative':float(max(abs(d))),'relative_range':[float(min(d)),float(max(d))],
        'max_difference_over_peak':float(max(abs(a-b))/max(b)),'worst_nm':float(lam[k])}
out={'ct_um':{'hard':3.75,'enge':3.75},'ds_m':1.875e-6,'frame':'lab','cases':{},
    'hard_tracking_policy':'Force RK4 ct0.9375um; exact retained nodes sampled at3.75um for radiation. Input sampling is checked against1.875/0.9375um.',
    'adaptive_policy':'Native adaptive is a cross-check; tighter Enge1km trials exceed55s. Production is fixed-step and independently halved.'}
for case,oldct in [('hard',.9375),('enge',1.875)]:
    for obs in ('11m','1km'):
        fine=get(case,obs,3.75,.9375e-6);main=get(case,obs,3.75,1.875e-6)
        mat=loadmat(O/f'MATLAB_pilot_{case}_{obs}_ct1.875um.mat',simplify_cells=True)['result']
        matxy=mat['W_radiation_lab'][:2].sum(0);matfull=mat['W_radiation_lab'].sum(0);matvel=mat['W_total_field_lab'][:2].sum(0)
        key=case+'_'+obs
        out['cases'][key]={'integration_vs_ds0p9375':{f'{p:g}':metrics(get(case,obs,3.75,p),fine) for p in (7.5e-6,3.75e-6,1.875e-6)},
            'input_sampling3p75_vs_finer':metrics(fine,get(case,obs,oldct,.9375e-6)),
            'native_vs_MATLAB_xy':metrics(main,matxy),'native_vs_MATLAB_total':metrics(main,matfull),
            'native_vs_LW_total_field_xy':metrics(main,matvel),'native_vs_analytic':metrics(main,ref[f'{obs}_S_analytic'][ids]),
            'MATLAB_z_fraction_max':float(max(mat['W_radiation_lab'][2]/matfull)),
            'velocity_inclusion_effect':metrics(matvel,matxy)}
        if case=='hard':
            out['cases'][key]['field_edge0p25_vs0p0625']=metrics(get(case,obs,.9375,.9375e-6,.25),get(case,obs,.9375,.9375e-6,.0625))
        auto=get(case,obs,3.75,1e-6,method=2)
        if auto is not None:out['cases'][key]['adaptive1e6_vs_manual_fine']=metrics(auto,fine)
        data=np.column_stack([lam,main,fine,matxy,matfull,matvel,ref[f'{obs}_S_analytic'][ids]])
        np.savetxt(O/f'pilot_{key}.txt',data,fmt='%.17e',header='Columns: lambda_nm SRW_ds1p875 SRW_ds0p9375 MATLAB_radiation_xy MATLAB_radiation_total MATLAB_radiation_plus_velocity_xy analytic')
        assert np.array_equal(np.loadtxt(O/f'pilot_{key}.txt'),data)
(O/'pilot_review.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
