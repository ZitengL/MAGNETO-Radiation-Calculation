"""Independent native adaptive-integrator spot checks, no production data replacement."""
import sys,json,time
sys.dont_write_bytecode=True
import calculate as m
import numpy as np
from srwpy.srwlib import SRWLWfr,SRWLPartBeam,srwl
precision=float(sys.argv[1]);tr=m.trajectory(1.875)
lam=np.load(m.O/'lambda_pilot.npy')[[0,10,20,30,50,70,89]]
fields=[];t=time.perf_counter()
for l in lam:
    w=SRWLWfr();w.allocate(1,1,1);w.mesh.eStart=w.mesh.eFin=m.HC/l
    w.mesh.xStart=w.mesh.xFin=w.mesh.yStart=w.mesh.yFin=0.;w.mesh.zStart=11
    w.partBeam=SRWLPartBeam(_Iavg=m.NATIVE_Q,_nPart=1,_partStatMom1=tr.partInitCond);w.unitElFld=1
    srwl.CalcElecFieldSR(w,tr,0,[2,precision,0,0,tr.np,1,0])
    fields.append([complex(w.arEx[0],w.arEx[1]),complex(w.arEy[0],w.arEy[1])])
f=np.asarray(fields).T
np.savez_compressed(m.O/f'adaptive_precision{precision:g}.npz',lambda_nm=lam,field_native=f,energy_components=abs(f)**2*m.HBAR*1e9*121,seconds=time.perf_counter()-t)
print(json.dumps({'precision':precision,'seconds':time.perf_counter()-t,'energy_xy':(abs(f)**2*m.HBAR*1e9*121).sum(axis=0).tolist()}),flush=True)
