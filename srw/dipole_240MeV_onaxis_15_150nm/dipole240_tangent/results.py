"""Assemble direct checkpoints only; component-resolved and analytic comparisons."""
import sys,json
sys.dont_write_bytecode=True
import model as m
import numpy as np
from scipy.io import loadmat
from scipy.special import kv
O=m.O
def collect(case,obs,grid='baseline',precision=1.875e-6,ct=3.75,edge=.125,b1=1e5):
    lam=np.load(O/f'lambda_{grid}.npy');label=f'{grid}_{case}_{obs}_lab_ct{ct:g}_edge{edge:g}_b{b1:g}_m0_p{precision:g}'
    folder=O/'blocks'/label;f=np.empty((2,len(lam)),complex);seen=np.zeros(len(lam),int);seconds=0.
    source=m.T/f'SRW_ct{ct:g}um_field2.5um.npz' if case=='enge' else O/f'hard_b{b1:g}_ct{ct:g}_edge{edge:g}.npz'
    sourcehash=m.track.sha(source);codehash=m.track.sha(O/'model.py');manifesthash=m.track.sha(O/'manifest.json')
    for path in sorted(folder.glob('*.npz')):
        if '.pending.' in path.name:continue
        a=np.load(path);c=json.loads(str(a['configuration']));i,j=c['start'],c['end']
        # Early10-node diagnostics coexist with full84-node pilots; prefer full pilots.
        if grid=='pilot' and (folder/f'000000_{len(lam):06d}.npz').exists() and j-i!=len(lam):continue
        assert np.array_equal(a['lambda_nm'],lam[i:j]);assert c['trajectory_sha256']==sourcehash
        assert c['code_sha256']==codehash and c['manifest_sha256']==manifesthash
        f[:,i:j]=a['field_native'];seen[i:j]+=1;seconds+=float(a['seconds'])
    assert np.all(seen==1),f'Incomplete or overlapping blocks: {folder}'
    assert np.all(np.isfinite(f))
    r=json.loads((O/'manifest.json').read_text())['observers'][obs]['R_source_m']
    return lam,f,abs(f)**2*m.HBAR*1e9*r**2,seconds

def metrics(a,b,lam):
    d=a/b-1;k=int(np.argmax(abs(d)))
    return {'max_abs_relative':float(max(abs(d))),'relative_range':[float(min(d)),float(max(d))],
        'max_difference_over_reference_peak':float(max(abs(a-b))/max(b)),
        'worst_wavelength_nm':float(lam[k]),'relative_at_first':float(d[0]),'relative_at_last':float(d[-1])}

def analytic(lam):
    c=json.loads((O/'manifest.json').read_text());q=1.6021892e-19;eps=8.85418781e-12
    omega=2*np.pi*m.C/(lam*1e-9);v=omega/c['omega_c']
    return 3*q*q*c['gamma_GPT']**2/(16*np.pi**3*eps*m.C)*v*v*kv(2/3,v/2)**2

def assemble():
    ref=np.load(O/'reference_original.npz');cfg=json.loads((O/'pilot_review.json').read_text())
    allsummary={'cases':{},'production_ds_m':cfg['ds_m'],'frame':'lab','ct_um':cfg['ct_um'],'bands_nm':[[15,25],[25,40],[40,70],[70,150]]}
    fieldchecks={}
    for name,b1 in [('hard',1e5),('hard_b1e6',1e6)]:
        g=np.load(O/f'GPT_{name}.npz')
        expected=m.analytical_field(g['x_m'],g['y_m'],g['z_m']-1,b1)
        actual=np.column_stack([g[k] for k in ('Bx_T','By_T','Bz_T')])
        fieldchecks[name]=float(np.max(abs(expected-actual)))
        assert fieldchecks[name]<1e-12,'Reconstructed field differs from actual GPT field'
    allsummary['analytic_GPT_field_formula_checks_max_abs_T']=fieldchecks
    manifest=json.loads((O/'manifest.json').read_text());totalCPU=0.
    for case in ('hard','enge'):
        for obs in ('11m','1km'):
            ct=cfg['ct_um'][case];ds=cfg['ds_m'];key=case+'_'+obs
            lam,f,pol,t1=collect(case,obs,'baseline',ds,ct);ld,fd,pd,t2=collect(case,obs,'dense',ds,ct)
            _,fc,pc,t3=collect(case,obs,'baseline',ds*2,ct);W=pol.sum(0);Wd=pd.sum(0);Wc=pc.sum(0);totalCPU+=t1+t2+t3
            assert np.all(pol[0]==0),'Broken x=0 orbit/observer polarization symmetry'
            assert np.array_equal(lam,ref['lambda_nm']) and np.array_equal(ld[::2],lam)
            main=loadmat(O/f'MATLAB_baseline_{case}_{obs}_ct1.875um.mat',simplify_cells=True)['result']
            original=loadmat(O/f'MATLAB_baseline_{case}_{obs}_ct3.75um.mat',simplify_cells=True)['result']
            dense=loadmat(O/f'MATLAB_dense_{case}_{obs}_ct1.875um.mat',simplify_cells=True)['result']
            assert np.array_equal(main['lambda_nm'],lam) and np.array_equal(dense['lambda_nm'],ld)
            old=ref[f'{obs}_S_MAGNETO_'+('hardedge' if case=='hard' else 'enge')]
            fineOriginal=ref[f'{obs}_S_'+('hard' if case=='hard' else 'enge')+'_dz1p875']
            ar=main['W_radiation_lab'];av=main['W_velocity_lab'];at=main['W_total_field_lab'];ai=main['W_interference_lab']
            xy=ar[:2].sum(0);tot=ar.sum(0);velxy=at[:2].sum(0);ad=dense['W_radiation_lab'];xd=ad[:2].sum(0)
            an=ref[f'{obs}_S_analytic'];ancheck=analytic(lam)
            summary={'native_vs_original_MAGNETO_total':metrics(W,old,lam),'native_vs_MAGNETO_fine_xy':metrics(W,xy,lam),
                'native_vs_LW_total_field_xy':metrics(W,velxy,lam),'native_vs_analytic_total':metrics(W,an,lam),
                'MAGNETO_original_vs_analytic':metrics(old,an,lam),'integral_main_vs_coarser':metrics(W,Wc,lam),
                'MATLAB_radiation_recomputed3p75_vs_original':metrics(original['W_radiation_lab'].sum(0),old,lam),
                'MATLAB_radiation_recomputed1p875_vs_provided_fine':metrics(tot,fineOriginal,lam),
                'MATLAB_reference_sampling3p75_vs1p875':metrics(original['W_radiation_lab'].sum(0),tot,lam),
                'MATLAB_velocity_inclusion_effect_xy':metrics(velxy,xy,lam),
                'MATLAB_z_fraction_range':[float(min(ar[2]/tot)),float(max(ar[2]/tot))],
                'analytic_formula_crosscheck_max_relative':float(max(abs(ancheck/an-1))),
                'dense_shared_complex_field_bitwise_equal':bool(np.array_equal(f,fd[:,::2])),
                'dense_shared_max_relative':float(max(abs(W/Wd[::2]-1))),
                'native_peak':float(max(W)),'native_peak_lambda_nm':float(lam[np.argmax(W)]),
                'native_x_fraction_of_peak':float(max(pol[0])/max(W)),
                'R_source_m':manifest['observers'][obs]['R_source_m'],'native_integrator_seconds':[t1,t2,t3]}
            midpoint=np.interp(np.log(ld[1::2]),np.log(lam),W)
            summary['wavelength_grid']={'max_midpoint_difference_over_peak':float(max(abs(Wd[1::2]-midpoint))/max(Wd)),
                'max_midpoint_pointwise_relative':float(max(abs(midpoint/Wd[1::2]-1)))}
            bands=[]
            for lo,hi in allsummary['bands_nm']:
                use=(lam>=lo-1e-12)&(lam<=hi+1e-12)
                bands.append({'range_nm':[lo,hi],'native_vs_original_total':metrics(W[use],old[use],lam[use]),
                    'native_vs_same_xy':metrics(W[use],xy[use],lam[use]),'native_vs_analytic':metrics(W[use],an[use],lam[use]),
                    'integration':metrics(W[use],Wc[use],lam[use])})
            summary['bands']=bands
            integrals=[]
            for lo,hi in [[15,150],*allsummary['bands_nm']]:
                row=[lo,hi]
                for ll,ww in [(lam,W),(ld,Wd),(lam,xy),(ld,xd)]:
                    use=(ll>lo)&(ll<hi);l=np.r_[lo,ll[use],hi];w=np.interp(np.log(l),np.log(ll),ww)
                    om=2*np.pi*m.C/(l*1e-9);row.append(float(-np.trapezoid(w,om)))
                row.extend([row[2]/row[3]-1,row[4]/row[5]-1]);integrals.append(row)
            integralarray=np.asarray(integrals)
            summary['band_energy_grid_max_relative_SRW']=float(max(abs(integralarray[:,-2])))
            summary['band_energy_grid_max_relative_MATLAB']=float(max(abs(integralarray[:,-1])))
            np.savetxt(O/f'band_grid_{key}.txt',integralarray,fmt='%.17e',header='Diagnostic integration over omega, band edges interpolated only for this check.\nColumns: lambda_lo_nm lambda_hi_nm SRW_baseline_J_sr SRW_dense_J_sr MATLAB_xy_baseline_J_sr MATLAB_xy_dense_J_sr SRW_baseline_div_dense_minus1 MATLAB_baseline_div_dense_minus1')
            assert np.array_equal(np.loadtxt(O/f'band_grid_{key}.txt'),integralarray)
            ids=np.unique(np.r_[np.linspace(0,4999,25).round().astype(int),np.arange(0,20,2),np.arange(4980,5000,2),
                np.argsort(abs(W/Wc-1))[-20:],np.argsort(abs(W/xy-1))[-10:]])
            grid='refinement_'+key
            np.save(O/f'lambda_{grid}.npy',lam[ids]);np.save(O/f'indices_{grid}.npy',ids);summary['refinement_nodes']=len(ids)
            np.savez_compressed(O/f'assembled_{key}.npz',lambda_nm=lam,omega=ref['omega'],field_native=f,energy_xy=W,energy_components=pol,
                energy_coarser=Wc,matlab_original_total=old,matlab_rad_xyz=ar,matlab_velocity_xyz=av,matlab_totalfield_xyz=at,
                matlab_interference_xyz=ai,analytic_total=an,dense_lambda_nm=ld,dense_field_native=fd,dense_energy_xy=Wd,
                dense_matlab_rad_xyz=ad,dense_matlab_totalfield_xyz=dense['W_total_field_lab'],dense_analytic_total=analytic(ld))
            allsummary['cases'][key]=summary
    allsummary['native_integrator_seconds_total']=totalCPU
    allsummary['production_wall_seconds']=json.loads((O/'full_completed.json').read_text())['wall_seconds']
    (O/'summary.json').write_text(json.dumps(allsummary,indent=2))
    print(json.dumps({k:{v:w for v,w in a.items() if v in ('native_vs_original_MAGNETO_total','native_vs_MAGNETO_fine_xy','integral_main_vs_coarser','wavelength_grid','refinement_nodes')} for k,a in allsummary['cases'].items()},indent=2),flush=True)

if __name__=='__main__':assemble()
