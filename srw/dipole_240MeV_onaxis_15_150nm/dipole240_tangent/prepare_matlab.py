"""Export original GPT arrays for new component-resolved MATLAB reference only."""
import sys
sys.dont_write_bytecode=True
import model as m
import numpy as np
from scipy.io import savemat
for case,step,name in [('hard',3.75,'GPT_hard.npz'),('hard',1.875,'GPT_hard_fine.npz'),
    ('enge',3.75,'GPT_enge_ct3.75um.npz'),('enge',1.875,'GPT_enge_ct1.875um.npz')]:
    a=np.load(m.O/name)
    tr={'t':a['time_s'],'r':np.stack([a[k] for k in ('x_m','y_m','z_m')]),
        'b':np.stack([a[k] for k in ('beta_x','beta_y','beta_z')]),'G':a['gamma'],
        'fB':np.stack([a[k] for k in ('Bx_T','By_T','Bz_T')]),'fE':np.zeros((3,len(a['ct_m']))),
        'q':-1.6021892e-19,'m':9.10953e-31,'nmacro':1.}
    savemat(m.O/f'GPT_{case}_ct{step:g}um.mat',{'tr':tr})
print('Original GPT arrays exported; no SRW trajectories used in MATLAB reference.')
