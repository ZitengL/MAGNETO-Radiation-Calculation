"""Native SRW six spectra at exact 240 MeV MATLAB spherical observer points."""
from pathlib import Path
from array import array
import sys,os,json,time,argparse
sys.dont_write_bytecode=True
O=Path(__file__).resolve().parent;BASE=O.parent;T=BASE/'undulator_240MeV_tracking'
sys.path.insert(0,str(BASE));import track
import numpy as np
from scipy.io import loadmat
from scipy.signal import find_peaks
from srwpy.srwlib import SRWLParticle,SRWLPrtTrj,SRWLWfr,SRWLPartBeam,srwl
os.environ.setdefault('MPLCONFIGDIR',str(O/'.matplotlib'))
C=299792458.;HBAR=1.05457182e-34;R=11.;NATIVE_Q=1.602189246e-19;NATIVE_HC=1000/.80654658
REF=track.SOURCE/'matlab/radiation_und15p_240MeV_single_electron/results/und240_spectra.mat'

def prepare():
    O.mkdir(exist_ok=True)
    trcfg=json.loads((T/'manifest.json').read_text());m=loadmat(REF,simplify_cells=True)['sp']
    lam=np.asarray(m['lambda_nm']).ravel();theta=np.asarray(m['theta_x_rad']).ravel()
    assert len(lam)==5000 and len(theta)==6 and np.max(abs(theta-np.arange(6)*1e-4))<1e-18
    pol=np.asarray(m['d2W_dwdOmega_xyz']);assert pol.shape==(6,3,5000)
    np.savez_compressed(O/'reference.npz',lambda_nm=lam,theta_rad=theta,omega=m['omega'],
        matlab_total=m['d2W_dwdOmega'],matlab_xy=pol[:,:2].sum(axis=1),matlab_z=pol[:,2],matlab_components=pol,
        matlab_convergence=m['conv_d2W_dwdOmega'],matlab_steps_um=m['conv_dz_um'],harmonics_nm=m['harmonics_nm'])
    dense=np.empty(2*len(lam)-1);dense[::2]=lam;dense[1::2]=np.sqrt(lam[:-1]*lam[1:])
    np.save(O/'lambda_baseline.npy',lam);np.save(O/'lambda_dense.npy',dense)
    ids=list(np.linspace(0,len(lam)-1,96).round().astype(int))
    for s in m['d2W_dwdOmega']:
        pk,_=find_peaks(s);pk=pk[np.argsort(s[pk])[-9:]]
        ids.extend(np.clip((pk[:,None]+np.arange(-2,3)).ravel(),0,len(lam)-1))
    ids=np.unique(ids);np.save(O/'pilot_indices.npy',ids);np.save(O/'lambda_pilot.npy',lam[ids])
    paths=[REF,Path(__file__),Path(track.__file__),T/'manifest.json',T/'tracking_summary.json',T/'field_rotated_compensated.npz']
    paths +=[T/f'SRW_ct{s:g}um.npz' for s in (60,30,15,7.5)]
    cfg={'scope':'240MeV undulator, six spectra only, 15..150nm',
        'gamma_GPT':trcfg['gamma_GPT'],'gamma_input':trcfg['gamma_SRW_interface'],'B_multiplier':trcfg['B_multiplier_all_components'],
        'theta_rad':theta.tolist(),'observer_rule':'D=11*(sin(theta),0,cos(theta)) m, exact per-angle spherical points',
        'baseline_nodes':len(lam),'dense_nodes':len(dense),'pilot_nodes_per_angle':len(ids),
        'normalization':'Native Iavg=1.602189246e-19 A; native alpha=1/137.0360411; SI |E_native|^2*hbar*1e9*11^2; hbar=1.05457182e-34',
        'phase':'native e_input_eV=(1000/0.80654658)/lambda_nm; no fitted wavelength shift',
        'comparison':'Native SRW Ex+Ey against original MATLAB x+y and total=x+y+z separately; no Ez reconstruction',
        'trajectory':'Reuse approved independent compensated trajectories. Full ct=0..1.369980000000408 m; no clipping, rotation change, extra drift or double compensation.',
        'tracking_manifest':str(T/'manifest.json'),'terminating':1,'radiation_method':0,
        'protected_SHA256':{str(p):track.sha(p) for p in paths},'tracking_input_SHA256':trcfg['SHA256']}
    (O/'manifest.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    print(json.dumps({'pilot_nodes':len(ids),'gamma_input':cfg['gamma_input'],'angles_mrad':(theta*1e3).tolist()}),flush=True)

def trajectory(step):
    a=np.load(T/f'SRW_ct{step:g}um.npz');cfg=json.loads((O/'manifest.json').read_text())
    assert float(a['gamma_input'])==cfg['gamma_input'] and float(a['field_multiplier'])==cfg['B_multiplier']
    tr=SRWLPrtTrj(_partInitCond=SRWLParticle(_gamma=float(a['gamma_input']),_relE0=1,_nq=-1))
    tr.allocate(len(a['ct_m']),True);tr.ctStart=float(a['ct_m'][0]);tr.ctEnd=float(a['ct_m'][-1])
    for k,v in track.COMPONENTS.items():setattr(tr,v,array('d',a[k]))
    return tr

def evaluate(tr,lam,theta,ds):
    energies=NATIVE_HC/np.asarray(lam);order=np.argsort(energies);e=energies[order]
    assert 1<=len(e)<=2
    w=SRWLWfr();w.allocate(len(e),1,1);w.mesh.eStart=float(e[0]);w.mesh.eFin=float(e[-1])
    w.mesh.xStart=w.mesh.xFin=float(R*np.sin(theta));w.mesh.yStart=w.mesh.yFin=0.;w.mesh.zStart=float(R*np.cos(theta))
    w.partBeam=SRWLPartBeam(_Iavg=NATIVE_Q,_nPart=1,_partStatMom1=tr.partInitCond);w.unitElFld=1
    srwl.CalcElecFieldSR(w,tr,0,[0,ds,0,0,tr.np,1,0])
    x=np.asarray(w.arEx,float).reshape(-1,2);y=np.asarray(w.arEy,float).reshape(-1,2)
    f=np.stack([x[:,0]+1j*x[:,1],y[:,0]+1j*y[:,1]])[:,np.argsort(order)]
    assert np.all(np.isfinite(f));return f

def run(grid,angle,step,ds,start,end):
    cfg=json.loads((O/'manifest.json').read_text());lam=np.load(O/f'lambda_{grid}.npy');end=min(end,len(lam))
    label=f'{grid}_a{angle}_ct{step:g}_ds{ds:g}';path=O/'blocks'/label/f'{start:06d}_{end:06d}.npz';path.parent.mkdir(parents=True,exist_ok=True)
    config={'grid':grid,'angle_index':angle,'step_um':step,'radiation_step_m':ds,'start':start,'end':end,
        'code_sha256':track.sha(__file__),'manifest_sha256':track.sha(O/'manifest.json'),
        'trajectory_sha256':track.sha(T/f'SRW_ct{step:g}um.npz')}
    if path.exists():assert json.loads(str(np.load(path)['configuration']))==config;return
    tr=trajectory(step);f=np.empty((2,end-start),complex);t0=time.perf_counter()
    for j in range(start,end,2):f[:,j-start:min(j+2,end)-start]=evaluate(tr,lam[j:min(j+2,end)],cfg['theta_rad'][angle],ds)
    elapsed=time.perf_counter()-t0;tmp=path.with_suffix('.pending.npz')
    np.savez_compressed(tmp,lambda_nm=lam[start:end],field_native=f,energy_components=abs(f)**2*HBAR*1e9*R**2,
        seconds=elapsed,configuration=json.dumps(config,sort_keys=True))
    tmp.replace(path);print(json.dumps({'label':label,'start':start,'end':end,'seconds':elapsed}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);p.add_argument('--grid',default='pilot')
    p.add_argument('--angle',type=int,default=0);p.add_argument('--step',type=float,default=7.5);p.add_argument('--ds',type=float,default=30e-6)
    p.add_argument('--start',type=int,default=0);p.add_argument('--end',type=int,default=100000);a=p.parse_args()
    if a.action=='prepare':prepare()
    else:run(a.grid,a.angle,a.step,a.ds,a.start,a.end)
