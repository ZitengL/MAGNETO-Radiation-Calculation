"""Complete dependent postprocessing once the current production job finishes."""
from pathlib import Path
import sys,time,json,subprocess
O=Path(__file__).resolve().parent
while not (O/'full_completed.json').exists():time.sleep(3)
full=json.loads((O/'full_completed.json').read_text())
assert full['failed']==0,'Production contains failed blocks'
guard=O.parents[1]/'srw_dipole_simulation/guarded_run.py'
for script,seconds,memory in [('analysis.py',60,768),('refine.py',None,None),('export.py',120,768),('write_readme.py',30,384)]:
    command=[sys.executable,'-B']
    if seconds is not None:command +=[str(guard),'--seconds',str(seconds),'--memory-mb',str(memory)]
    command +=[str(O/script)]
    with (O/f'pipeline_{Path(script).stem}.log').open('w',encoding='utf-8') as log:
        result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,creationflags=0x08000000)
    status={'stage':script,'returncode':result.returncode,'time_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    (O/'pipeline_status.json').write_text(json.dumps(status));print(json.dumps(status),flush=True)
    if result.returncode:raise SystemExit(result.returncode)
print('Production, refinement, figures, TXT validation and README complete; final visual/MATLAB review remains.',flush=True)
