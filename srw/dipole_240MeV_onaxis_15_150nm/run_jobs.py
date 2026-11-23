"""Four independent bounded SRW workers; checkpoints retain raw complex fields."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import sys,subprocess,json,time,argparse
O=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('stage',choices=['pilot','full']);a=p.parse_args();jobs=[]
if a.stage=='pilot':
    import calculate as m
    n=len(m.np.load(O/'lambda_pilot.npy'))
    settings=[(1.875,2.5,ds,'full',1) for ds in (3.75e-6,1.875e-6,.9375e-6,.46875e-6)]
    settings +=[(ct,2.5,.46875e-6,'full',1) for ct in (3.75,.9375)]
    settings +=[(.9375,1.25,.46875e-6,'full',1),(1.875,2.5,.46875e-6,'gpt',1),(1.875,2.5,.46875e-6,'full',0)]
    for ct,fs,ds,win,term in settings:jobs.append(('pilot',ct,fs,ds,0,n,win,term))
else:
    cfg=json.loads((O/'pilot_review.json').read_text());ct=cfg['production_ct_um'];fs=cfg['production_field_um'];ds=cfg['production_ds_m']
    for grid,n,d in [('baseline',5000,ds),('baseline',5000,ds*2),('dense',9999,ds)]:
        for start in range(0,n,64):jobs.append((grid,ct,fs,d,start,min(start+64,n),'full',1))

def work(job):
    grid,ct,fs,ds,start,end,win,term=job
    log=O/'logs'/f'{grid}_ct{ct:g}_field{fs:g}_ds{ds:g}_{win}_t{term}_{start:06d}.log';log.parent.mkdir(exist_ok=True)
    cmd=[sys.executable,'-B',str(O.parents[1]/'srw_dipole_simulation/guarded_run.py'),'--seconds','55','--memory-mb','640',str(O/'calculate.py'),
        'run','--grid',grid,'--ct',str(ct),'--field',str(fs),'--ds',str(ds),'--start',str(start),'--end',str(end),'--window',win,'--term',str(term)]
    t0=time.perf_counter()
    with log.open('w',encoding='utf-8') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,creationflags=0x08000000)
    return {'job':job,'returncode':r.returncode,'wall_seconds':time.perf_counter()-t0,'log':str(log)}
t0=time.perf_counter();done=[]
with ThreadPoolExecutor(max_workers=4) as pool:
    for f in as_completed([pool.submit(work,j) for j in jobs]):
        done.append(f.result());status={'completed':len(done),'total':len(jobs),'failed':sum(x['returncode']!=0 for x in done),'wall_seconds':time.perf_counter()-t0}
        (O/f'{a.stage}_progress.json').write_text(json.dumps(status));print(json.dumps(done[-1]),flush=True)
result={'workers':4,'wall_seconds':time.perf_counter()-t0,'failed':sum(x['returncode']!=0 for x in done),'jobs':done}
(O/f'{a.stage}_completed.json').write_text(json.dumps(result,indent=2))
if result['failed']:raise SystemExit(1)
