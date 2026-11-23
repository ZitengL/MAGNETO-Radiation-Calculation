"""Checks selected after inspecting all original-grid integration and benchmark differences."""
from concurrent.futures import ThreadPoolExecutor,as_completed
import sys,json,time,subprocess
from pathlib import Path
O=Path(__file__).resolve().parent
import calculate as m
n=len(m.np.load(O/'lambda_refinement.npy'))
settings=[(1.875,2.5,.1171875e-6),(.9375,2.5,.234375e-6),(.9375,1.25,.234375e-6)]
jobs=[(ct,fs,ds,i,min(i+20,n)) for ct,fs,ds in settings for i in range(0,n,20)]
def work(job):
    ct,fs,ds,i,j=job;log=O/'logs'/f'refine_ct{ct:g}_field{fs:g}_ds{ds:g}_{i:06d}.log'
    cmd=[sys.executable,'-B',str(O.parents[1]/'srw_dipole_simulation/guarded_run.py'),'--seconds','55','--memory-mb','768',str(O/'calculate.py'),
        'run','--grid','refinement','--ct',str(ct),'--field',str(fs),'--ds',str(ds),'--start',str(i),'--end',str(j)]
    with log.open('w',encoding='utf-8') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,creationflags=0x08000000)
    return {'job':job,'returncode':r.returncode,'log':str(log)}
t=time.perf_counter();done=[]
with ThreadPoolExecutor(max_workers=4) as pool:
    for f in as_completed([pool.submit(work,j) for j in jobs]):done.append(f.result());print(json.dumps(done[-1]),flush=True)
result={'wall_seconds':time.perf_counter()-t,'failed':sum(x['returncode']!=0 for x in done),'jobs':done}
(O/'refinement_completed.json').write_text(json.dumps(result,indent=2))
if result['failed']:raise SystemExit(1)
