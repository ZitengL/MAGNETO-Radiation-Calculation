from pathlib import Path
import sys,time,json,subprocess
O=Path(__file__).resolve().parent
while not (O/'full_completed.json').exists():time.sleep(3)
assert json.loads((O/'full_completed.json').read_text())['failed']==0
guard=O.parents[2]/'srw_dipole_simulation/guarded_run.py'
for script,seconds,memory in [('results.py',90,768),('refine.py',None,None),('export.py',150,1024),('write_readme.py',30,384)]:
    cmd=[sys.executable,'-B']
    if seconds is not None:cmd +=[str(guard),'--seconds',str(seconds),'--memory-mb',str(memory)]
    cmd +=[str(O/script)]
    with (O/f'pipeline_{Path(script).stem}.log').open('w',encoding='utf-8') as log:
        result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,creationflags=0x08000000)
    status={'stage':script,'returncode':result.returncode,'UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    (O/'pipeline_status.json').write_text(json.dumps(status));print(json.dumps(status),flush=True)
    if result.returncode:raise SystemExit(result.returncode)
print('Tangent four-case computation/export complete; final image/MATLAB review remains.',flush=True)
