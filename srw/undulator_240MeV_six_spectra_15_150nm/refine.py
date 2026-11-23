"""Check worst full-grid integration discrepancies with finer native steps; keep production unchanged."""
import analyze as a
import calculate as m
import numpy as np,json,subprocess,sys,time
from concurrent.futures import ThreadPoolExecutor,as_completed
O=m.O
jobs=[];selected={}
for i in range(6):
    b=a.collect('baseline',i,7.5,30e-6);c=a.collect('baseline',i,7.5,60e-6);l=b['lambda_nm']
    ref=np.load(O/'reference.npz');from scipy.ndimage import maximum_filter1d
    half=round(.01/np.mean(np.diff(np.log(l))));env=maximum_filter1d(b['energy'],size=2*half+1)
    error=abs(b['energy']-c['energy'])/env;ids=[]
    for lo,hi in [(15,30),(30,60),(60,110),(110,150.0001)]:
        use=np.flatnonzero((l>=lo*(1-1e-12))&(l<hi));ids.extend(use[np.argsort(error[use])[-12:]])
    for h in ref['harmonics_nm'][i]:
        k=int(np.argmin(abs(l-h)));ids.extend(range(max(0,k-2),min(len(l),k+3)))
    ids=np.unique(np.r_[ids,0,len(l)-1]);selected[i]=ids
    grid=f'check_angle{i}';np.save(O/f'lambda_{grid}.npy',l[ids]);np.save(O/f'check_indices_angle{i}.npy',ids)
    for ct,ds in [(7.5,30e-6),(7.5,15e-6),(7.5,7.5e-6),(15,30e-6)]:
        jobs.append((grid,i,ct,ds,0,len(ids)))
def work(job):
    grid,i,ct,ds,start,end=job;log=O/'logs'/f'{grid}_ct{ct:g}_ds{ds:g}.log'
    cmd=[sys.executable,'-B',str(O.parents[1]/'srw_dipole_simulation/guarded_run.py'),'--seconds','55','--memory-mb','640',str(O/'calculate.py'),
        'run','--grid',grid,'--angle',str(i),'--step',str(ct),'--ds',str(ds),'--end',str(end)]
    with log.open('w',encoding='utf-8') as stream:r=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,creationflags=0x08000000)
    assert r.returncode==0,str(log);return job
t0=time.perf_counter();done=[]
with ThreadPoolExecutor(max_workers=2) as pool:
    for f in as_completed([pool.submit(work,j) for j in jobs]):done.append(f.result())
summary={}
for i,ids in selected.items():
    grid=f'check_angle{i}';base=a.collect('baseline',i,7.5,30e-6);coarse=a.collect('baseline',i,7.5,60e-6)
    cols=[base['lambda_nm'][ids],base['energy'][ids],coarse['energy'][ids]];row={}
    for ct,ds in [(7.5,30e-6),(7.5,15e-6),(7.5,7.5e-6),(15,30e-6)]:
        q=a.collect(grid,i,ct,ds);cols.append(q['energy'])
        row[f'ct{ct:g}_ds{ds:g}']={'max_diff_over_primary_peak':float(max(abs(q['energy']-base['energy'][ids]))/max(base['energy']))}
    np.savetxt(O/f'refined_checks_theta_{i/10:.1f}mrad.txt',np.column_stack(cols),fmt='%.17e',header='Selected from full-grid integration discrepancies in four wavelength bands and original harmonic markers; no formal values replaced.\n'
        'Columns: lambda_nm primary_native_ct7p5_rad30um control_native_ct7p5_rad60um repeated_primary native_rad15um native_rad7p5um native_ct15_rad30um\nEnergy units J/(rad/s)/sr; no fitted shifts or normalization.')
    summary[str(i)]={'nodes':len(ids),'comparisons':row}
(O/'refinement_summary.json').write_text(json.dumps({'angles':summary,'jobs':len(done),'wall_seconds':time.perf_counter()-t0},indent=2))
print(json.dumps(summary,indent=2),flush=True)
