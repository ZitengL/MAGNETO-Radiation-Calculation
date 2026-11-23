"""Isolate finer FORCE tracking at fixed3.75um radiation input sampling."""
from array import array
import sys,json,time
sys.dont_write_bytecode=True
import model as m
import numpy as np
observer=sys.argv[1]
source=m.O/'hard_b100000_ct0.46875_edge0.125.npz';a=np.load(source);stride=8;ct=a['ct_m'][::stride]
cfg=json.loads((m.O/'manifest.json').read_text());p=m.SRWLParticle(_gamma=cfg['gamma_input'],_relE0=1,_nq=-1)
tr=m.SRWLPrtTrj(_partInitCond=p);tr.allocate(len(ct),True);tr.ctStart=float(ct[0]);tr.ctEnd=float(ct[-1])
for key,attr in m.track.COMPONENTS.items():setattr(tr,attr,array('d',a[key][::stride]))
lam=np.load(m.O/'lambda_pilot.npy');f=np.empty((2,len(lam)),complex);t=time.perf_counter()
for j in range(0,len(lam),2):f[:,j:j+2]=m.evaluate(tr,lam[j:j+2],observer,'lab',0,.9375e-6)
r=cfg['observers'][observer]['R_source_m']
np.savez_compressed(m.O/f'force_refined_{observer}.npz',lambda_nm=lam,field_native=f,energy_components=abs(f)**2*m.HBAR*1e9*r*r,
    force_ct_um=.46875,radiation_input_ct_um=3.75,ds_m=.9375e-6,source_sha256=m.track.sha(source),seconds=time.perf_counter()-t)
print(json.dumps({'observer':observer,'force_ct_um':.46875,'input_ct_um':3.75,'seconds':time.perf_counter()-t}),flush=True)
