"""Resumable direct native SRW calls with four memory/time-bounded workers."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import sys,subprocess,json,time,argparse
O=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('stage',choices=['pilot','full']);a=p.parse_args();jobs=[]
if a.stage=='pilot':
    import calculate as m
    n=len(m.np.load(O/'lambda_pilot.npy'))
    settings=[(angle,7.5,ds) for angle in range(6) for ds in (60e-6,30e-6,15e-6)]
    settings +=[(angle,step,15e-6) for angle in (0,3,5) for step in (60,30,15)]
    for angle,step,ds in settings:
        for start in range(0,n,96):jobs.append(('pilot',angle,step,ds,start,min(start+96,n)))
else:
    cfg=json.loads((O/'pilot_review.json').read_text());ds=cfg['production_radiation_step_m'];ct=cfg['production_ct_step_um']
    for grid,n,d in [('baseline',5000,ds),('baseline',5000,ds*2),('dense',9999,ds)]:
        for angle in range(6):
            for start in range(0,n,256):jobs.append((grid,angle,ct,d,start,min(start+256,n)))

def work(job):
    grid,angle,ct,ds,start,end=job;log=O/'logs'/f'{grid}_a{angle}_ct{ct:g}_ds{ds:g}_{start:06d}.log';log.parent.mkdir(exist_ok=True)
    cmd=[sys.executable,'-B',str(O.parents[1]/'srw_dipole_simulation/guarded_run.py'),'--seconds','55','--memory-mb','640',str(O/'calculate.py'),
        'run','--grid',grid,'--angle',str(angle),'--step',str(ct),'--ds',str(ds),'--start',str(start),'--end',str(end)]
    t0=time.perf_counter()
    with log.open('w',encoding='utf-8') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,creationflags=0x08000000)
    return {'job':job,'returncode':r.returncode,'wall_seconds':time.perf_counter()-t0,'log':str(log)}
t0=time.perf_counter();done=[]
with ThreadPoolExecutor(max_workers=4) as pool:
    for f in as_completed([pool.submit(work,j) for j in jobs]):
        done.append(f.result());status={'completed':len(done),'total':len(jobs),'failed':sum(x['returncode']!=0 for x in done),'wall_seconds':time.perf_counter()-t0}
        try:(O/f'{a.stage}_progress.json').write_text(json.dumps(status))
        except PermissionError:pass
        print(json.dumps(done[-1]),flush=True)
out={'workers':4,'wall_seconds':time.perf_counter()-t0,'failed':sum(x['returncode']!=0 for x in done),'jobs':done}
(O/f'{a.stage}_completed.json').write_text(json.dumps(out,indent=2))
if out['failed']:raise SystemExit(1)
