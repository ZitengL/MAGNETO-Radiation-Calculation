"""Four bounded native SRW workers with atomic checkpoints."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import sys,json,time,subprocess,argparse
O=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('stage',choices=['pilot','sampling','full']);args=p.parse_args();jobs=[]
if args.stage=='pilot':
    import model as m
    n=len(m.np.load(O/'lambda_pilot.npy'))
    for case,ct in [('hard',.9375),('enge',1.875)]:
        for observer in ('11m','1km'):
            for ds in (7.5e-6,3.75e-6,1.875e-6,.9375e-6):jobs.append(('pilot',case,observer,ct,.125,1e5,0,ds,0,n))
            for prec in (1e-5,1e-6):jobs.append(('pilot',case,observer,ct,.125,1e5,2,prec,0,n))
    for observer in ('11m','1km'):
        jobs +=[('pilot','hard',observer,.46875,.125,1e5,0,.9375e-6,0,n),
            ('pilot','hard',observer,.9375,.25,1e5,0,.9375e-6,0,n),
            ('pilot','hard',observer,.9375,.0625,1e5,0,.9375e-6,0,n),
            ('pilot','enge',observer,.9375,.125,1e5,0,.9375e-6,0,n)]
elif args.stage=='sampling':
    import model as m
    n=len(m.np.load(O/'lambda_pilot.npy'))
    for case in ('hard','enge'):
        for obs in ('11m','1km'):
            for ds in (7.5e-6,3.75e-6,1.875e-6,.9375e-6):jobs.append(('pilot',case,obs,3.75,.125,1e5,0,ds,0,n))
            jobs.append(('pilot',case,obs,3.75,.125,1e5,2,1e-6,0,n))
else:
    cfg=json.loads((O/'pilot_review.json').read_text())
    for case in ('hard','enge'):
        ct=cfg['ct_um'][case]
        for obs in ('11m','1km'):
            for grid,n,ds in [('baseline',5000,cfg['ds_m']),('baseline',5000,cfg['ds_m']*2),('dense',9999,cfg['ds_m'])]:
                for i in range(0,n,128):jobs.append((grid,case,obs,ct,.125,1e5,0,ds,i,min(i+128,n)))
def work(job):
    grid,case,obs,ct,edge,b1,method,prec,i,j=job
    log=O/'logs'/f'{grid}_{case}_{obs}_ct{ct:g}_edge{edge:g}_m{method}_p{prec:g}_{i:06d}.log';log.parent.mkdir(exist_ok=True)
    cmd=[sys.executable,'-B',str(O.parents[2]/'srw_dipole_simulation/guarded_run.py'),'--seconds','55','--memory-mb','1024',str(O/'model.py'),
        'run','--grid',grid,'--case',case,'--observer',obs,'--ct',str(ct),'--edge',str(edge),'--b1',str(b1),'--method',str(method),'--precision',str(prec),'--start',str(i),'--end',str(j)]
    t=time.perf_counter()
    with log.open('w',encoding='utf-8') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,creationflags=0x08000000)
    return {'job':job,'returncode':r.returncode,'wall_seconds':time.perf_counter()-t,'log':str(log)}
t=time.perf_counter();done=[]
with ThreadPoolExecutor(max_workers=4) as pool:
    for f in as_completed([pool.submit(work,j) for j in jobs]):
        done.append(f.result());status={'completed':len(done),'total':len(jobs),'failed':sum(v['returncode']!=0 for v in done),'wall_seconds':time.perf_counter()-t}
        (O/f'{args.stage}_progress.json').write_text(json.dumps(status));print(json.dumps(done[-1]),flush=True)
result={'workers':4,'failed':sum(v['returncode']!=0 for v in done),'wall_seconds':time.perf_counter()-t,'jobs':done}
(O/f'{args.stage}_completed.json').write_text(json.dumps(result,indent=2))
if result['failed']:raise SystemExit(1)
