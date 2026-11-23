"""Pinned GPT/MATLAB tangent benchmark; new outputs only, native SRW binary unchanged."""
from pathlib import Path
from array import array
import sys,os,json,time,argparse,importlib.util
sys.dont_write_bytecode=True
O=Path(__file__).resolve().parent;BASE=O.parents[1];T=BASE/'dipole_240MeV_tracking'
sys.path.insert(0,str(BASE));import track
import numpy as np
from scipy.io import loadmat,savemat
from srwpy.srwlib import SRWLParticle,SRWLPrtTrj,SRWLMagFldC,SRWLWfr,SRWLPartBeam,srwl
os.environ.setdefault('MPLCONFIGDIR',str(O/'.matplotlib'))
SRC=track.SOURCE;REF=SRC/'matlab/dipole_analytic_comparison/dipole240_tangent'
HARD=SRC/'gpt_sims/single_electron/E240MeV_dipole_L10cm_hardedge'
B_SCALE=(1.6021892e-19/9.10953e-31)/(1.602176462e-19/9.10938188e-31)
C=299792458.;HBAR=1.05457182e-34;NATIVE_Q=1.602189246e-19;HC=1000/.80654658
_spec=importlib.util.spec_from_file_location('existing_reader',BASE/'undulator_240MeV_tracking/tracking.py')
_reader=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_reader)

def prepare():
    O.mkdir(exist_ok=True)
    d=loadmat(REF/'results/compare_dipole_tangent.mat',simplify_cells=True)['out']
    lam=np.asarray(d['lambda_nm']).ravel();assert len(lam)==5000
    P=np.asarray(d['tangent_point']).ravel();direc=np.asarray(d['tangent_dir']).ravel()
    theta=float(d['bend_rad'])/2;cs=np.cos(theta);sn=np.sin(theta)
    rot=np.array([[1,0,0],[0,cs,-sn],[0,sn,cs]])
    ref={'lambda_nm':lam,'omega':d['omega'],'tangent_point':P,'tangent_direction':direc,'rotation_lab_to_tangent':rot}
    observers={}
    for name,key in [('11m','o11m'),('1km','o1km')]:
        a=d[key];obs=np.asarray(a['obs']).ravel();r=float(a['R_src'])
        assert abs(np.linalg.norm(obs-P)/r-1)<1e-14
        observers[name]={'position_lab_m':obs.tolist(),'position_tangent_m':(rot@(obs-P)).tolist(),'R_source_m':r,
            'R_origin_m':float(np.linalg.norm(obs)),'original_MAGNETO_normalization_factor':float((r/np.linalg.norm(obs))**2)}
        for field in ('S_analytic','S_MAGNETO_hardedge','S_MAGNETO_enge','S_hard_dz1p875','S_hard_b1e6','S_enge_dz1p875'):
            ref[f'{name}_{field}']=np.asarray(a[field]).ravel()
    np.savez_compressed(O/'reference_original.npz',**ref)
    dense=np.empty(9999);dense[::2]=lam;dense[1::2]=np.sqrt(lam[:-1]*lam[1:])
    ids=np.unique(np.r_[np.linspace(0,4999,65).round().astype(int),np.arange(0,20,2),np.arange(4980,5000,2)])
    for name,v in [('baseline',lam),('dense',dense),('pilot',lam[ids])]:np.save(O/f'lambda_{name}.npy',v)
    np.save(O/'pilot_indices.npy',ids)
    trajectories={}
    paths=[REF/'README.md',REF/'compare_dipole_tangent.m',REF/'results/compare_dipole_tangent.mat',
        HARD/'README.md',HARD/'single_electron_dipole_hardedge.in',
        Path(os.environ.get('GPT_ROOT', 'GPT'), r'elems\sources\BStatic\rectmagnet.c'),track.SRW/'env/python/srwpy/srwlpy.pyd',
        T/'manifest.json',T/'SRW_ct1.875um_field2.5um.npz',T/'SRW_ct0.9375um_field1.25um.npz',Path(__file__),Path(track.__file__)]
    for name,path in [('hard',HARD/'result.gdf'),('hard_fine',HARD/'tout_convergence/result_dz1p875um.gdf'),
        ('hard_b1e6',HARD/'b1_check/result_b1e6.gdf')]:
        a,p=_reader.gpt_reference(path);np.savez_compressed(O/f'GPT_{name}.npz',**a)
        trajectories[name]={'path':str(path),'samples':len(a['ct_m']),'ct_range_m':[float(a['ct_m'][0]),float(a['ct_m'][-1])],
            'gamma':float(a['gamma'][0]),**p};paths.append(path)
    old=json.loads((T/'manifest.json').read_text());oldgdf=old['reference']['3.75']['path']
    assert track.sha(oldgdf)==old['SHA256'][oldgdf],'Original Enge GPT trajectory changed; re-evaluate before reusing tracking'
    for step in (3.75,1.875):
        a=np.load(T/f'GPT_ct{step:g}um.npz');np.savez_compressed(O/f'GPT_enge_ct{step:g}um.npz',**{k:a[k] for k in a.files})
    gamma=float(d['gamma']);assert abs(gamma-old['gamma_GPT'])<1e-10
    cfg={'scope':'User confirmed both hard/Enge models and both 11m/1km observers',
        'gamma_GPT':gamma,'gamma_input':gamma*.511003414/.51099890221,'B_multiplier':B_SCALE,
        'B0_T':float(d['B0']),'rho_m':float(d['rho']),'omega_c':float(d['omega_c']),
        'bend_rad':float(d['bend_rad']),'tangent_theta_rad':theta,'P_lab_m':P.tolist(),'direction_lab':direc.tolist(),
        'rotation_lab_to_tangent':rot.tolist(),'observers':observers,'GPT':trajectories,
        'Enge_approved_end_ct_m':old['reference']['3.75']['ct_range_m'][1],
        'normalization':'|E_native|^2*hbar*1e9*R_source^2; original MATLAB norm(origin)^2 is already converted to R_source^2 in target MAT; do not apply distance factor twice',
        'native_binary_unchanged':True,'no_fitted_correction':'No fitted spectral shifts, amplitude corrections, smoothing, or component reconstruction.',
        'protected_SHA256':{str(p):track.sha(p) for p in paths},'original_Enge_tracking_SHA256':old['SHA256']}
    (O/'manifest.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    savemat(O/'reference_geometry.mat',{'lambda_nm':lam,'lambda_pilot_nm':lam[ids],'P':P,'rotation':rot,
        'obs11m':np.asarray(observers['11m']['position_lab_m']),'obs1km':np.asarray(observers['1km']['position_lab_m'])})
    print(json.dumps({k:v for k,v in cfg.items() if 'SHA256' not in k},indent=2),flush=True)

def analytical_field(x,y,z,b1):
    x,y,z=np.broadcast_arrays(x,y,z);d=abs(z)-.05;e=np.exp(b1*d);h=-b1*x
    bn=.65/(1+2*e*np.cos(h)+e*e)
    b=np.stack([-bn*(1+e*np.cos(h)),np.zeros_like(x),-bn*e*np.sin(h)*np.where(z<0,-1.,1.)],axis=-1)
    inside=(abs(y)<=.1)&(d*b1<=10)
    assert not np.any(inside&(abs(x)*b1>=np.pi)&(d<=0))
    b[~inside]=0;return b

def hard_field(b1,edge_um):
    # Three adjacent uniform tables. All points use the full GPT formula.
    # This avoids assuming a different physical plateau model or using nonuniform native lookup.
    outer=.05+10/b1;inner=.05-40/b1
    n=round((outer-inner)/(edge_um*1e-6))
    zlist=[np.linspace(-outer,-inner,n+1),np.linspace(-inner,inner,3),np.linspace(inner,outer,n+1)]
    x=np.array([-.1/b1,0,.1/b1]);y=np.array([-.1,0,.1]);fields=[];saved={'x':x,'y':y,'b1':b1,'edge_requested_um':edge_um,'B_multiplier':B_SCALE}
    for j,z in enumerate(zlist):
        b=analytical_field(x[None,None,:],y[None,:,None],z[:,None,None],b1)*B_SCALE
        one=track.field_container(x,y,z,b);fields.append(one.arMagFld[0]);saved[f'z{j}']=z;saved[f'B{j}']=b
    mag=SRWLMagFldC(fields,array('d',[0]*3),array('d',[0]*3),array('d',[1]*3))
    return mag,saved

def track_hard(ctstep,edge,b1):
    cfg=json.loads((O/'manifest.json').read_text());end=cfg['GPT']['hard']['ct_range_m'][1]
    mag,f=hard_field(b1,edge)
    p=SRWLParticle(_gamma=cfg['gamma_input'],_relE0=1,_nq=-1);n=round(end/(ctstep*1e-6))+1
    tr=SRWLPrtTrj(_partInitCond=p);tr.allocate(n,True);tr.ctStart=0;tr.ctEnd=end
    t=time.perf_counter();srwl.CalcPartTraj(tr,mag,[1]);seconds=time.perf_counter()-t
    a={'ct_m':np.linspace(0,end,n),**{k:np.asarray(getattr(tr,v)).copy() for k,v in track.COMPONENTS.items()}}
    assert all(np.all(np.isfinite(v)) for v in a.values()) and np.all(a['x_m']==0) and np.all(a['beta_x']==0)
    label=f'hard_b{b1:g}_ct{ctstep:g}_edge{edge:g}'
    np.savez_compressed(O/(label+'.npz'),**a,gamma_input=p.gamma,gamma_target=cfg['gamma_GPT'],field_multiplier=B_SCALE)
    np.savez_compressed(O/f'field_b{b1:g}_edge{edge:g}.npz',**f)
    g=np.load(O/('GPT_hard.npz' if b1==1e5 else 'GPT_hard_b1e6.npz'));c=g['ct_m']
    info={'ct_step_um':ctstep,'edge_step_um':edge,'b1':b1,'points':n,'tracking_seconds':seconds,
        'final_bend_mrad':float(np.arctan2(a['beta_y'][-1],a['beta_z'][-1])*1e3),
        'vs_GPT_at_actual_ct':{k:float(max(abs(np.interp(c,a['ct_m'],a[k])-g[k]))) for k in list(track.COMPONENTS)[:6]}}
    (O/(label+'.json')).write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)

def trajectory(case='enge',frame='lab',step=1.875,edge=.125,b1=1e5):
    cfg=json.loads((O/'manifest.json').read_text())
    path=T/f'SRW_ct{step:g}um_field2.5um.npz' if case=='enge' else O/f'hard_b{b1:g}_ct{step:g}_edge{edge:g}.npz'
    a=np.load(path);assert float(a['gamma_input'])==cfg['gamma_input'] and float(a['field_multiplier'])==B_SCALE
    p=SRWLParticle(_gamma=cfg['gamma_input'],_relE0=1,_nq=-1)
    arrays={k:a[k] for k in track.COMPONENTS};ct=a['ct_m']
    if frame=='tangent':
        # Exact rigid coordinate transformation and change of time origin at a stored node.
        # This is a diagnostic of native representation dependence, never a fitted trajectory shift.
        R=np.asarray(cfg['rotation_lab_to_tangent']);P=np.asarray(cfg['P_lab_m'])
        pos=R@(np.stack([arrays[k] for k in ('x_m','y_m','z_m')])-P[:,None])
        vel=R@np.stack([arrays[k] for k in ('beta_x','beta_y','beta_z')])
        bf=R@np.stack([arrays[k] for k in ('Bx_T','By_T','Bz_T')])
        j=int(np.argmin(abs(a['z_m']-P[2])))
        ct=ct-ct[j]
        for keys,mat in [(('x_m','y_m','z_m'),pos),(('beta_x','beta_y','beta_z'),vel),(('Bx_T','By_T','Bz_T'),bf)]:
            arrays.update({k:mat[i] for i,k in enumerate(keys)})
        p.x,p.y,p.z=map(float,pos[:,j]);p.xp,p.yp=map(float,vel[:2,j])
    tr=SRWLPrtTrj(_partInitCond=p);tr.allocate(len(ct),True);tr.ctStart=float(ct[0]);tr.ctEnd=float(ct[-1])
    for k,v in track.COMPONENTS.items():setattr(tr,v,array('d',arrays[k]))
    return tr,path

def evaluate(tr,lam,observer,frame,method,precision,term=1):
    cfg=json.loads((O/'manifest.json').read_text());obs=cfg['observers'][observer]['position_'+('lab' if frame=='lab' else 'tangent')+'_m']
    energy=HC/np.asarray(lam);order=np.argsort(energy);e=energy[order]
    w=SRWLWfr();w.allocate(len(e),1,1);w.mesh.eStart=float(e[0]);w.mesh.eFin=float(e[-1])
    w.mesh.xStart=w.mesh.xFin=float(obs[0]);w.mesh.yStart=w.mesh.yFin=float(obs[1]);w.mesh.zStart=float(obs[2])
    w.partBeam=SRWLPartBeam(_Iavg=NATIVE_Q,_nPart=1,_partStatMom1=tr.partInitCond);w.unitElFld=1
    srwl.CalcElecFieldSR(w,tr,0,[method,precision,0,0,tr.np,term,0])
    ex=np.asarray(w.arEx,float).reshape(-1,2);ey=np.asarray(w.arEy,float).reshape(-1,2)
    f=np.stack([ex[:,0]+1j*ex[:,1],ey[:,0]+1j*ey[:,1]])[:,np.argsort(order)]
    assert np.all(np.isfinite(f));return f

def run(case,observer,frame,grid,step,edge,b1,method,precision,start,end):
    cfg=json.loads((O/'manifest.json').read_text());lam=np.load(O/f'lambda_{grid}.npy');end=min(end,len(lam))
    label=f'{grid}_{case}_{observer}_{frame}_ct{step:g}_edge{edge:g}_b{b1:g}_m{method}_p{precision:g}'
    dest=O/'blocks'/label/f'{start:06d}_{end:06d}.npz';dest.parent.mkdir(parents=True,exist_ok=True)
    tr,path=trajectory(case,frame,step,edge,b1)
    conf={'case':case,'observer':observer,'frame':frame,'grid':grid,'ct_um':step,'edge_um':edge,'b1':b1,
        'method':method,'precision':precision,'start':start,'end':end,'code_sha256':track.sha(__file__),
        'manifest_sha256':track.sha(O/'manifest.json'),'trajectory_sha256':track.sha(path)}
    if dest.exists():assert json.loads(str(np.load(dest)['configuration']))==conf;return
    t=time.perf_counter();f=np.empty((2,end-start),complex)
    for j in range(start,end,2):f[:,j-start:min(j+2,end)-start]=evaluate(tr,lam[j:min(j+2,end)],observer,frame,method,precision)
    seconds=time.perf_counter()-t;radius=cfg['observers'][observer]['R_source_m'];tmp=dest.with_suffix('.pending.npz')
    np.savez_compressed(tmp,lambda_nm=lam[start:end],field_native=f,energy_components=abs(f)**2*HBAR*1e9*radius**2,
        configuration=json.dumps(conf,sort_keys=True),seconds=seconds);tmp.replace(dest)
    print(json.dumps({'label':label,'start':start,'end':end,'seconds':seconds}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','track','run']);p.add_argument('--case',default='enge',choices=['enge','hard'])
    p.add_argument('--observer',default='11m',choices=['11m','1km']);p.add_argument('--frame',default='lab',choices=['lab','tangent']);p.add_argument('--grid',default='pilot')
    p.add_argument('--ct',type=float,default=1.875);p.add_argument('--edge',type=float,default=.125);p.add_argument('--b1',type=float,default=1e5)
    p.add_argument('--method',type=int,default=0);p.add_argument('--precision',type=float,default=1.875e-6);p.add_argument('--start',type=int,default=0);p.add_argument('--end',type=int,default=100000)
    a=p.parse_args()
    if a.action=='prepare':prepare()
    elif a.action=='track':track_hard(a.ct,a.edge,a.b1)
    else:run(a.case,a.observer,a.frame,a.grid,a.ct,a.edge,a.b1,a.method,a.precision,a.start,a.end)
