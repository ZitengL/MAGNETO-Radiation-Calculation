"""Test exact binary increments in the unchanged native manual integrator.
Only an exactly-zero transverse incoming integrand is omitted, with full trajectory supplied.
All endpoints chosen from trajectory metadata before evaluating any spectrum, never fitted.
"""
import sys,time,json
sys.dont_write_bytecode=True
import calculate as m
import numpy as np
from srwpy.srwlib import SRWLWfr,SRWLPartBeam,srwl
power=int(sys.argv[1]);ds=2.**(-power);tr=m.trajectory(1.875)
native_end=(tr.ctEnd-tr.ctStart)/(tr.np-1)*(tr.np-1)
base=2.**-20;intervals=2*int(np.ceil((native_end-.87)/base/2));start=native_end-intervals*base
n=intervals*2**(power-20)
assert (native_end-start)/ds==n and n%2==0 and (native_end-start)/n==ds
g=np.load(m.T/'SRW_ct1.875um_field2.5um.npz');mask=g['ct_m']<=start
assert all(np.all(g[k][mask]==0) for k in ('x_m','y_m','beta_x','beta_y','Bx_T','By_T','Bz_T'))
# Since start and native_end have magnitude <2m and share a 2**-52 lattice,
# adding ds>=2**-24 is exact over this interval (no repeated-addition drift).
lam=np.load(m.O/'lambda_pilot.npy');ids=np.array([0,2,6,10,20,30,50,70,89]);lam=lam[ids]
fields=[];t=time.perf_counter()
for l in lam:
    w=SRWLWfr();w.allocate(1,1,1);w.mesh.eStart=w.mesh.eFin=m.HC/l
    w.mesh.xStart=w.mesh.xFin=w.mesh.yStart=w.mesh.yFin=0.;w.mesh.zStart=11
    w.partBeam=SRWLPartBeam(_Iavg=m.NATIVE_Q,_nPart=1,_partStatMom1=tr.partInitCond);w.unitElFld=1
    srwl.CalcElecFieldSR(w,tr,0,[0,ds,start,native_end,tr.np,1,0])
    fields.append([complex(w.arEx[0],w.arEx[1]),complex(w.arEy[0],w.arEy[1])])
f=np.asarray(fields).T;energy=abs(f)**2*m.HBAR*1e9*121
np.savez_compressed(m.O/f'binary_step_2powminus{power}.npz',lambda_nm=lam,indices=ids,field_native=f,energy_components=energy,
    seconds=time.perf_counter()-t,ds=ds,start=start,end=native_end,intervals=n)
print(json.dumps({'power':power,'seconds':time.perf_counter()-t,'start':start,'end':native_end,'energy_xy':energy.sum(0).tolist()}),flush=True)
