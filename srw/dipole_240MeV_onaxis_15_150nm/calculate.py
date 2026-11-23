"""Direct native SRW radiation from the approved independent 240 MeV dipole tracking."""
from pathlib import Path
from array import array
import sys,os,json,time,argparse
sys.dont_write_bytecode=True
O=Path(__file__).resolve().parent;BASE=O.parent;T=BASE/'dipole_240MeV_tracking'
sys.path.insert(0,str(BASE));import track
import numpy as np
from scipy.io import loadmat,savemat
from srwpy.srwlib import SRWLParticle,SRWLPrtTrj,SRWLWfr,SRWLPartBeam,srwl
C=299792458.;HBAR=1.05457182e-34;R=11.;NATIVE_Q=1.602189246e-19;HC=1000/.80654658
REF=track.SOURCE/'matlab/radiation_dipole240_single_electron/results/dipole240_onaxis.mat'

def prepare():
    O.mkdir(exist_ok=True)
    tc=json.loads((T/'manifest.json').read_text());d=loadmat(REF,simplify_cells=True)['d']
    lam=np.asarray(d['lambda_nm']).ravel();pol=np.asarray(d['d2W_dwdOmega_xyz'])
    assert len(lam)==5000 and pol.shape==(3,5000)
    np.savez_compressed(O/'reference_original.npz',lambda_nm=lam,omega=d['omega'],energy_xyz=pol,
        energy_xy=pol[:2].sum(axis=0),energy_total=d['d2W_dwdOmega'],photon=d['d2N_dwdOmega'],
        conv_energy=d['conv_d2W_dwdOmega'],conv_steps_um=d['conv_dz_um'])
    dense=np.empty(9999);dense[::2]=lam;dense[1::2]=np.sqrt(lam[:-1]*lam[1:])
    ids=np.unique(np.r_[np.round(np.linspace(0,4999,81)).astype(int),np.arange(0,20,2)])
    for name,a in [('baseline',lam),('dense',dense),('pilot',lam[ids])]:np.save(O/f'lambda_{name}.npy',a)
    np.save(O/'pilot_indices.npy',ids)
    for step in (3.75,1.875):
        g=np.load(T/f'GPT_ct{step:g}um.npz')
        tr={'t':g['time_s'],'r':np.stack([g[k] for k in ('x_m','y_m','z_m')]),
            'b':np.stack([g[k] for k in ('beta_x','beta_y','beta_z')]),'G':g['gamma'],
            'fB':np.stack([g[k] for k in ('Bx_T','By_T','Bz_T')]),'fE':np.zeros((3,len(g['ct_m']))),
            'q':tc['reference'][f'{step:g}']['q'],'m':tc['reference'][f'{step:g}']['m'],'nmacro':1.}
        savemat(O/f'GPT_ct{step:g}um.mat',{'tr':tr})
    paths=[REF,Path(__file__),Path(track.__file__),T/'manifest.json',T/'tracking_summary.json',
        T/'field_2.5um.npz',T/'field_1.25um.npz',track.SRW/'env/python/srwpy/srwlpy.pyd',
        track.SOURCE/'magneto/emission_time_spectrum.m']
    paths += list(T.glob('SRW_*.npz'))+list(T.glob('GPT_*.npz'))
    cfg={'scope':'240 MeV short dipole, on axis at (0,0,11)m, 15..150nm, 5000 original nodes + 9999 independently evaluated dense nodes',
        'gamma_GPT':tc['gamma_GPT'],'gamma_input':tc['gamma_SRW_interface'],'B_multiplier':tc['B_multiplier'],
        'native_photon_input':'(1000/0.80654658)/lambda_nm; exact native wavelength convention, no fitted shift',
        'normalization':'native Iavg=1.602189246e-19 A; |E_native|^2*hbar*1e9*11^2, hbar=1.05457182e-34; native alpha unchanged',
        'comparison':'native Ex+Ey vs original MATLAB x+y and original x+y+z; new exact LW acceleration+velocity reference, complex fields added before squaring; no Ez reconstruction',
        'trajectory':'approved compensated SRW tracking; full original ct=0..1.1499975000025333m, no added drift, no repeated B or gamma compensation',
        'native_model_limit':'Native SRW reconstructs longitudinal phase using its paraxial trajectory representation; supplied arZ is not an exact-geometry phase override. This remains a native implementation difference at 81mrad, not an added external approximation.',
        'pilot_nodes':len(ids),'protected_SHA256':{str(p):track.sha(p) for p in paths},'tracking_original_SHA256':tc['SHA256']}
    (O/'manifest.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    print(json.dumps({'pilot_nodes':len(ids),'gamma_input':cfg['gamma_input']}),flush=True)

def trajectory(step,fs=2.5,window='full'):
    a=np.load(T/f'SRW_ct{step:g}um_field{fs:g}um.npz');cfg=json.loads((O/'manifest.json').read_text())
    assert float(a['gamma_input'])==cfg['gamma_input'] and float(a['field_multiplier'])==cfg['B_multiplier']
    sl=slice(None)
    if window=='gpt':
        start=json.loads((T/'manifest.json').read_text())['reference']['3.75']['ct_range_m'][0]
        first=int(np.argmin(abs(a['ct_m']-start)));sl=slice(first,None)
    # partInitCond remains the original ct=0 emission condition, also for a sub-window.
    # Native SRW uses ct+partInitCond.z as its longitudinal integration coordinate.
    p=SRWLParticle(_gamma=cfg['gamma_input'],_relE0=1,_nq=-1)
    tr=SRWLPrtTrj(_partInitCond=p);tr.allocate(len(a['ct_m'][sl]),True)
    tr.ctStart=float(a['ct_m'][sl][0]);tr.ctEnd=float(a['ct_m'][sl][-1])
    for k,v in track.COMPONENTS.items():setattr(tr,v,array('d',a[k][sl]))
    return tr

def evaluate(tr,lam,ds,term=1):
    order=np.argsort(HC/lam);e=(HC/lam)[order]
    w=SRWLWfr();w.allocate(len(e),1,1);w.mesh.eStart=float(e[0]);w.mesh.eFin=float(e[-1])
    w.mesh.xStart=w.mesh.xFin=w.mesh.yStart=w.mesh.yFin=0.;w.mesh.zStart=R
    w.partBeam=SRWLPartBeam(_Iavg=NATIVE_Q,_nPart=1,_partStatMom1=tr.partInitCond);w.unitElFld=1
    srwl.CalcElecFieldSR(w,tr,0,[0,ds,0,0,tr.np,term,0])
    x=np.asarray(w.arEx,float).reshape(-1,2);y=np.asarray(w.arEy,float).reshape(-1,2)
    f=np.stack([x[:,0]+1j*x[:,1],y[:,0]+1j*y[:,1]])[:,np.argsort(order)]
    assert np.all(np.isfinite(f));return f

def run(grid,ct,fs,ds,start,end,window='full',term=1):
    lam=np.load(O/f'lambda_{grid}.npy');end=min(end,len(lam))
    label=f'{grid}_ct{ct:g}_field{fs:g}_ds{ds:g}_{window}_t{term}'
    dest=O/'blocks'/label/f'{start:06d}_{end:06d}.npz';dest.parent.mkdir(parents=True,exist_ok=True)
    conf={'grid':grid,'ct_um':ct,'field_um':fs,'ds_m':ds,'start':start,'end':end,'window':window,'terminating':term,
        'code_sha256':track.sha(__file__),'manifest_sha256':track.sha(O/'manifest.json'),
        'trajectory_sha256':track.sha(T/f'SRW_ct{ct:g}um_field{fs:g}um.npz')}
    if dest.exists():assert json.loads(str(np.load(dest)['configuration']))==conf;return
    tr=trajectory(ct,fs,window);f=np.empty((2,end-start),complex);t0=time.perf_counter()
    for j in range(start,end,2):f[:,j-start:min(j+2,end)-start]=evaluate(tr,lam[j:min(j+2,end)],ds,term)
    sec=time.perf_counter()-t0;tmp=dest.with_suffix('.pending.npz')
    np.savez_compressed(tmp,lambda_nm=lam[start:end],field_native=f,energy_components=abs(f)**2*HBAR*1e9*121,
        seconds=sec,configuration=json.dumps(conf,sort_keys=True));tmp.replace(dest)
    print(json.dumps({'label':label,'start':start,'end':end,'seconds':sec}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);p.add_argument('--grid',default='pilot')
    p.add_argument('--ct',type=float,default=1.875);p.add_argument('--field',type=float,default=2.5);p.add_argument('--ds',type=float,default=.9375e-6)
    p.add_argument('--start',type=int,default=0);p.add_argument('--end',type=int,default=100000);p.add_argument('--window',choices=['full','gpt'],default='full');p.add_argument('--term',type=int,default=1)
    a=p.parse_args()
    if a.action=='prepare':prepare()
    else:run(a.grid,a.ct,a.field,a.ds,a.start,a.end,a.window,a.term)
