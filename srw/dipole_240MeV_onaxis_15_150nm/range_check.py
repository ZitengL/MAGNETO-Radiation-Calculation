"""Native integration start test: omit only an exactly zero transverse integrand.
The full ct=0 trajectory and its original initial condition remain supplied.
No trajectory reslicing, re-origining, endpoint extension, phase fit or field change.
"""
import sys,time,json
sys.dont_write_bytecode=True
import calculate as m
import numpy as np
from srwpy.srwlib import SRWLWfr,SRWLPartBeam,srwl
ds=float(sys.argv[1]);start=float(sys.argv[2]);tr=m.trajectory(1.875)
lam=np.load(m.O/'lambda_pilot.npy');ids=np.array([0,10,20,30,50,70,89]);lam=lam[ids]
# Directly verify the omitted transverse trajectory and field are exactly zero.
a=np.load(m.T/'SRW_ct1.875um_field2.5um.npz');mask=a['ct_m']<=start
assert all(np.all(a[k][mask]==0) for k in ('x_m','y_m','beta_x','beta_y','Bx_T','By_T','Bz_T'))
fields=[];t=time.perf_counter()
for l in lam:
    w=SRWLWfr();w.allocate(1,1,1);w.mesh.eStart=w.mesh.eFin=m.HC/l
    w.mesh.xStart=w.mesh.xFin=w.mesh.yStart=w.mesh.yFin=0.;w.mesh.zStart=11
    w.partBeam=SRWLPartBeam(_Iavg=m.NATIVE_Q,_nPart=1,_partStatMom1=tr.partInitCond);w.unitElFld=1
    srwl.CalcElecFieldSR(w,tr,0,[0,ds,start,tr.ctEnd,tr.np,1,0])
    fields.append([complex(w.arEx[0],w.arEx[1]),complex(w.arEy[0],w.arEy[1])])
f=np.asarray(fields).T
np.savez_compressed(m.O/f'range_explicit_end_start{start:g}_ds{ds:g}.npz',lambda_nm=lam,indices=ids,field_native=f,energy_components=abs(f)**2*m.HBAR*1e9*121,seconds=time.perf_counter()-t)
print(json.dumps({'start_m':start,'ds':ds,'seconds':time.perf_counter()-t,'energy_xy':(abs(f)**2*m.HBAR*1e9*121).sum(axis=0).tolist()}),flush=True)
