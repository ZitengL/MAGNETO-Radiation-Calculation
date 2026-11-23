"""Collect direct evaluations and measure native integration/trajectory stability."""
import calculate as m
import numpy as np,json
from scipy.ndimage import maximum_filter1d
O=m.O

def collect(grid,angle,ct,ds):
    lam=np.load(O/f'lambda_{grid}.npy');field=np.empty((2,len(lam)),complex);covered=np.zeros(len(lam),bool);seconds=0;count=0
    folder=O/'blocks'/f'{grid}_a{angle}_ct{ct:g}_ds{ds:g}'
    for p in sorted(folder.glob('*.npz')):
        a=np.load(p);cfg=json.loads(str(a['configuration']));lo=cfg['start'];hi=cfg['end']
        assert not np.any(covered[lo:hi]);assert np.array_equal(a['lambda_nm'],lam[lo:hi]);assert cfg['manifest_sha256']==m.track.sha(O/'manifest.json')
        field[:,lo:hi]=a['field_native'];covered[lo:hi]=True;seconds+=float(a['seconds']);count+=1
    assert np.all(covered),f'{folder}: missing {sum(~covered)}'
    energy_components=abs(field)**2*m.HBAR*1e9*m.R**2
    return {'lambda_nm':lam,'field':field,'components':energy_components,'energy':energy_components.sum(axis=0),'seconds':seconds,'blocks':count}

def measure(lam,a,b):
    # Same +/-1% logarithmic-wavelength envelope as original MATLAB convergence analysis.
    half=round(.01/float(np.mean(np.diff(np.log(lam)))))
    env=maximum_filter1d(np.asarray(b),size=2*half+1,mode='nearest')
    d=abs(a-b);nz=env>0
    return {'max_absolute_difference':float(max(d)),'max_difference_over_peak':float(max(d)/max(b)),
        'max_difference_over_local_envelope':float(max(d[nz]/env[nz])),
        'max_pointwise_relative_difference':float(max(d[b>0]/b[b>0])),'zero_reference_nodes':int(sum(b==0))}

def pilot():
    cfg=json.loads((O/'manifest.json').read_text());ref=np.load(O/'reference.npz');ids=np.load(O/'pilot_indices.npy')
    out={};timings=[]
    # Sparse sample spacings are irregular. Envelope from full MATLAB spectrum is used for pilot comparisons.
    lam=ref['lambda_nm'];half=round(.01/np.mean(np.diff(np.log(lam))))
    for angle in range(6):
        runs={d:collect('pilot',angle,7.5,d) for d in (60e-6,30e-6,15e-6)};fine=runs[15e-6]['energy'];env=maximum_filter1d(ref['matlab_xy'][angle],size=2*half+1)[ids]
        row={}
        for d,a in runs.items():
            row[str(d)]={'difference_vs15um_over_MATLAB_peak':float(max(abs(a['energy']-fine))/max(ref['matlab_xy'][angle])),
                'difference_vs15um_over_MATLAB_local_envelope':float(max(abs(a['energy']-fine)/env)),
                'seconds':a['seconds']};timings.append(a['seconds']/len(fine))
        row['native15um_vs_MATLAB_xy']={'max_difference_over_peak':float(max(abs(fine-ref['matlab_xy'][angle,ids]))/max(ref['matlab_xy'][angle])),
            'max_difference_over_local_envelope':float(max(abs(fine-ref['matlab_xy'][angle,ids])/env))}
        if angle in (0,3,5):
            row['trajectory_checks']={}
            for ct in (60,30,15):
                a=collect('pilot',angle,ct,15e-6)
                row['trajectory_checks'][str(ct)]={'difference_over_MATLAB_peak':float(max(abs(a['energy']-fine))/max(ref['matlab_xy'][angle])),
                    'difference_over_MATLAB_local_envelope':float(max(abs(a['energy']-fine)/env))}
        out[str(angle)]=row
    result={'pilot_results':out,'production_ct_step_um':7.5,'production_radiation_step_m':30e-6,
        'full_integration_control_step_m':60e-6,'finer_sample_step_m':15e-6,
        'mean_pilot_seconds_per_point':float(np.mean(timings)),
        'note':'Pilot envelope uses full original MATLAB grid. Full comparisons later use each native reference full-grid envelope. No accuracy inferred at zero crossings.'}
    (O/'pilot_review.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':pilot()
