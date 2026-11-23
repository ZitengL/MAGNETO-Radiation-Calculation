"""Preserve accurate RK4 integration while independently checking radiation input sampling."""
import sys,json
sys.dont_write_bytecode=True
import model as m
import numpy as np
source=m.O/'hard_b100000_ct0.9375_edge0.125.npz';a=np.load(source)
for target in (1.875,3.75):
    stride=round(target/.9375);data={k:a[k][::stride] for k in ['ct_m',*m.track.COMPONENTS]}
    assert len(data['ct_m'])==(len(a['ct_m'])-1)//stride+1 and data['ct_m'][-1]==a['ct_m'][-1]
    path=m.O/f'hard_b100000_ct{target:g}_edge0.125.npz'
    assert not path.exists(),'Do not overwrite an independently tracked trajectory'
    np.savez_compressed(path,**data,gamma_input=a['gamma_input'],gamma_target=a['gamma_target'],field_multiplier=a['field_multiplier'],
        integration_ct_um=.9375,saved_sampling_ct_um=target,decimation_stride=stride,source_sha256=m.track.sha(source))
    info={'purpose':'Radiation input sampling check only; exact retained nodes, no interpolation or retracking',
        'force_integrator_ct_um':.9375,'radiation_input_sampling_ct_um':target,'source':str(source),'source_sha256':m.track.sha(source),'stride':stride}
    path.with_suffix('.json').write_text(json.dumps(info,indent=2))
    print(json.dumps(info),flush=True)
