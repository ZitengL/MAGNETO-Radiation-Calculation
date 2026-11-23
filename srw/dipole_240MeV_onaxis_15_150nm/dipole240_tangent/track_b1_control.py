"""Native RK4 b1=1e6 control; exact original free drift before a local tracking origin.
This rebases time inside an already specified field-free region, not an external trajectory extension.
"""
import sys,time,json
sys.dont_write_bytecode=True
import model as m
import numpy as np
cfg=json.loads((m.O/'manifest.json').read_text());b1=1e6;edge=.0125;h=.09375e-6;origin_ct=.93
end=cfg['GPT']['hard']['ct_range_m'][1];beta0=np.sqrt(1-cfg['gamma_GPT']**-2)
mag,f=m.hard_field(b1,edge)
p=m.SRWLParticle(_gamma=cfg['gamma_input'],_relE0=1,_nq=-1,_z=beta0*origin_ct)
n=round((end-origin_ct)/h)+1;tr=m.SRWLPrtTrj(_partInitCond=p);tr.allocate(n,True);tr.ctStart=0;tr.ctEnd=end-origin_ct
t=time.perf_counter();m.srwl.CalcPartTraj(tr,mag,[1]);seconds=time.perf_counter()-t
localct=np.linspace(origin_ct,end,n);ct=np.linspace(0,end,round(end/3.75e-6)+1);post=ct>=origin_ct
a={'ct_m':ct}
for key,attr in m.track.COMPONENTS.items():
    initial=beta0*ct if key=='z_m' else np.full(len(ct),beta0) if key=='beta_z' else np.zeros(len(ct))
    initial[post]=np.interp(ct[post],localct,np.asarray(getattr(tr,attr)));a[key]=initial
assert np.all(a['x_m']==0) and np.all(a['beta_x']==0)
label=f'hard_b{b1:g}_ct3.75_edge{edge:g}'
np.savez_compressed(m.O/(label+'.npz'),**a,gamma_input=cfg['gamma_input'],gamma_target=cfg['gamma_GPT'],field_multiplier=m.B_SCALE,
    integration_ct_um=h*1e6,saved_sampling_ct_um=3.75,local_origin_global_ct_m=origin_ct)
np.savez_compressed(m.O/f'field_b{b1:g}_edge{edge:g}.npz',**f)
g=np.load(m.O/'GPT_hard_b1e6.npz');c=g['ct_m']
info={'b1':b1,'edge_um':edge,'force_ct_um':h*1e6,'radiation_sampling_ct_um':3.75,
    'time_origin_global_ct_m':origin_ct,'tracking_seconds':seconds,
    'exact_free_drift':'From original ct0 toct0.93, before the first nonzero field. r=(0,0,beta0*ct), beta=(0,0,beta0), no force or model extrapolation.',
    'resampling':'Linear evaluation of native local output onto original global uniform3.75um nodes; no fitted shift.',
    'final_bend_mrad':float(np.arctan2(a['beta_y'][-1],a['beta_z'][-1])*1e3),
    'vs_GPT_at_actual_ct':{k:float(max(abs(np.interp(c,ct,a[k])-g[k]))) for k in list(m.track.COMPONENTS)[:6]}}
(m.O/(label+'.json')).write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2),flush=True)
