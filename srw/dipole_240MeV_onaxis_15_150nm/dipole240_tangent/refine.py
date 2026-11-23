"""Selected finer-radiation and finer-trajectory checks after full-grid comparison."""
from concurrent.futures import ThreadPoolExecutor,as_completed
import sys,json,time,subprocess
sys.dont_write_bytecode=True
import model as m
O=m.O;summary=json.loads((O/'summary.json').read_text());jobs=[]
for case in ('hard','enge'):
    for obs in ('11m','1km'):
        grid='refinement_'+case+'_'+obs;n=summary['cases'][case+'_'+obs]['refinement_nodes']
        for ct,ds in [(3.75,.9375e-6),(1.875,1.875e-6)]:
            for i in range(0,n,24):jobs.append((grid,case,obs,ct,ds,i,min(i+24,n)))
def work(job):
    grid,case,obs,ct,ds,i,j=job;log=O/'logs'/f'{grid}_ct{ct:g}_p{ds:g}_{i:06d}.log'
    cmd=[sys.executable,'-B',str(O.parents[2]/'srw_dipole_simulation/guarded_run.py'),'--seconds','55','--memory-mb','1024',str(O/'model.py'),
        'run','--grid',grid,'--case',case,'--observer',obs,'--ct',str(ct),'--precision',str(ds),'--start',str(i),'--end',str(j)]
    t=time.perf_counter()
    with log.open('w',encoding='utf-8') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,creationflags=0x08000000)
    return {'job':job,'returncode':r.returncode,'wall_seconds':time.perf_counter()-t,'log':str(log)}
t=time.perf_counter();done=[]
with ThreadPoolExecutor(max_workers=4) as pool:
    for f in as_completed([pool.submit(work,j) for j in jobs]):done.append(f.result());print(json.dumps(done[-1]),flush=True)
result={'failed':sum(v['returncode']!=0 for v in done),'wall_seconds':time.perf_counter()-t,'jobs':done}
(O/'refinement_completed.json').write_text(json.dumps(result,indent=2))
if result['failed']:raise SystemExit(1)
