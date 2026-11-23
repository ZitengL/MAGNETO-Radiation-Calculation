"""Assemble only direct native outputs; compare unsmoothed spectra and independent grids."""
import sys,json
from pathlib import Path
sys.dont_write_bytecode=True
import calculate as m
import numpy as np
from scipy.io import loadmat
O=m.O
def collect(grid,ds,ct=1.875,fs=2.5):
    lam=np.load(O/f'lambda_{grid}.npy');folder=O/'blocks'/f'{grid}_ct{ct:g}_field{fs:g}_ds{ds:g}_full_t1'
    out=np.empty((2,len(lam)),complex);seen=np.zeros(len(lam),int);seconds=0.
    for p in sorted(folder.glob('*.npz')):
        if '.pending.' in p.name:continue
        a=np.load(p);c=json.loads(str(a['configuration']));i,j=c['start'],c['end']
        assert np.array_equal(a['lambda_nm'],lam[i:j]);assert c['code_sha256']==m.track.sha(O/'calculate.py')
        assert c['manifest_sha256']==m.track.sha(O/'manifest.json')
        assert c['trajectory_sha256']==m.track.sha(m.T/f'SRW_ct{ct:g}um_field{fs:g}um.npz')
        out[:,i:j]=a['field_native'];seen[i:j]+=1;seconds+=float(a['seconds'])
    assert np.all(seen==1),f'{folder}: missing/overlapping nodes'
    assert np.all(np.isfinite(out))
    return lam,out,seconds

def metrics(a,b,lam):
    d=a-b;rel=d/b;k=int(np.argmax(abs(rel)))
    return {'max_abs_relative':float(max(abs(rel))),'p95_abs_relative':float(np.quantile(abs(rel),.95)),
        'max_abs_difference_over_reference_peak':float(max(abs(d))/max(b)),
        'signed_relative_at_first_node':float(rel[0]),'signed_relative_at_last_node':float(rel[-1]),
        'worst_relative_wavelength_nm':float(lam[k]),'reference_fraction_of_peak_at_worst':float(b[k]/max(b))}

def main():
    cfg=json.loads((O/'pilot_review.json').read_text());ds=cfg['production_ds_m'];ct=cfg['production_ct_um'];fs=cfg['production_field_um']
    lam,f,cpu1=collect('baseline',ds,ct,fs);ld,fd,cpu2=collect('dense',ds,ct,fs);_,fc,cpu3=collect('baseline',ds*2,ct,fs)
    pol=abs(f)**2*m.HBAR*1e9*121;W=pol.sum(0);Wd=(abs(fd)**2*m.HBAR*1e9*121).sum(0);Wc=(abs(fc)**2*m.HBAR*1e9*121).sum(0)
    assert np.array_equal(ld[::2],lam)
    # Shared wavelengths are independently evaluated, using different blocks/pairings.
    same=bool(np.array_equal(f,fd[:,::2]));shared=float(max(abs(W-Wd[::2]))/max(W))
    ref=np.load(O/'reference_original.npz');assert np.array_equal(ref['lambda_nm'],lam)
    mat=loadmat(O/'matlab_baseline_ct1.875um.mat',simplify_cells=True)['result']
    mat3=loadmat(O/'matlab_baseline_ct3.75um.mat',simplify_cells=True)['result']
    md=loadmat(O/'matlab_dense_ct1.875um.mat',simplify_cells=True)['result']
    assert np.array_equal(mat['lambda_nm'],lam) and np.array_equal(md['lambda_nm'],ld)
    old=ref['energy_xy'];oldtotal=ref['energy_total'];new=mat['W_total_field'][:2].sum(0);new3=mat3['W_total_field'][:2].sum(0)
    newd=md['W_total_field'][:2].sum(0)
    radiation3=mat3['W_radiation'][:2].sum(0)
    cpu=cpu1+cpu2+cpu3
    summary={'status':'all original and dense wavelengths directly computed; preserve numerical sensitivity at weakest short wavelengths',
        'baseline_nodes':len(lam),'dense_nodes':len(ld),'production_ct_um':ct,'production_field_um':fs,'production_ds_m':ds,
        'native_vs_original_matlab_xy':metrics(W,old,lam),'native_vs_original_matlab_total':metrics(W,oldtotal,lam),
        'native_vs_velocity_inclusive_reference_xy':metrics(W,new,lam),
        'integration_ds_vs_twice_ds':metrics(W,Wc,lam),'MATLAB_total_field_ct3p75_vs1p875':metrics(new3,new,lam),
        'new_radiation_reference_vs_original_radiation':metrics(radiation3,old,lam),
        'dense_shared_complex_fields_bitwise_equal':same,'dense_shared_max_difference_over_peak':shared,
        'native_energy_density_peak':float(max(W)),'native_peak_lambda_nm':float(lam[np.argmax(W)]),
        'native_15nm_fraction_of_peak':float(W[0]/max(W)),'native_x_max_fraction':float(max(pol[0])/max(W)),
        'original_matlab_z_max_fraction':float(max(ref['energy_xyz'][2]/ref['energy_total'])),
        'new_LW_z_max_fraction':float(max(mat['W_total_field'][2]/mat['W_total_field'].sum(0))),
        'native_integrator_seconds_sum':cpu,'integrator_seconds_baseline_dense_control':[cpu1,cpu2,cpu3],
        'production_wall_seconds':json.loads((O/'full_completed.json').read_text())['wall_seconds'],
        'normalization_only_expected_relative_offset':-15.990660e-6,
        'no_rescaling':'All native SRW normalization retained; no fitted spectral shift, scaling, smoothing, sample replacement or discarded points.'}
    bands=[]
    for lo,hi in [(15,20),(20,30),(30,50),(50,100),(100,150)]:
        ix=(lam>=lo-1e-12)&(lam<=hi+1e-12)
        bands.append({'range_nm':[lo,hi],'native_vs_original':metrics(W[ix],old[ix],lam[ix]),
            'native_vs_velocity_inclusive':metrics(W[ix],new[ix],lam[ix]),'integration':metrics(W[ix],Wc[ix],lam[ix])})
    summary['bands']=bands
    def gridcheck(l,w,wfine):
        # This interpolation is ONLY a grid-resolution diagnostic, never production data.
        mid=np.interp(np.log(ld[1::2]),np.log(l),w)
        err=mid-wfine[1::2]
        return {'max_midpoint_linear_loglambda_reconstruction_error_over_peak':float(max(abs(err))/max(wfine)),
            'max_midpoint_pointwise_relative':float(max(abs(err/wfine[1::2])))}
    summary['wavelength_grid_SRW']=gridcheck(lam,W,Wd);summary['wavelength_grid_MATLAB_total_field']=gridcheck(lam,new,newd)
    ints=[]
    for lo,hi in [(15,150),(15,30),(30,60),(60,110),(110,150)]:
        row=[lo,hi]
        # Band integrals of dW/domega over omega, with exact band edges interpolated for diagnostics only.
        for l,w in [(lam,W),(ld,Wd),(lam,new),(ld,newd)]:
            inside=(l>lo)&(l<hi);ll=np.r_[lo,l[inside],hi];ww=np.interp(np.log(ll),np.log(l),w)
            om=2*np.pi*m.C/(ll*1e-9);row.append(float(-np.trapezoid(ww,om)))
        row +=[row[2]/row[3]-1,row[4]/row[5]-1];ints.append(row)
    summary['band_integral_max_relative_grid_change_SRW']=float(max(abs(np.asarray(ints)[:,-2])))
    summary['band_integral_max_relative_grid_change_MATLAB']=float(max(abs(np.asarray(ints)[:,-1])))
    np.savez_compressed(O/'assembled.npz',lambda_nm=lam,omega=ref['omega'],field_native=f,energy_components=pol,
        energy_xy=W,energy_coarser=Wc,dense_lambda_nm=ld,dense_field_native=fd,dense_energy_xy=Wd,
        matlab_original_xy=old,matlab_original_total=oldtotal,matlab_original_z=ref['energy_xyz'][2],
        matlab_total_field_xy=new,matlab_total_field_xyz=mat['W_total_field'],matlab_radiation_xyz=mat['W_radiation'],
        matlab_velocity_xyz=mat['W_velocity'],matlab_interference_xyz=mat['W_interference'],matlab_total_field_coarse=new3,
        dense_matlab_total_field_xy=newd)
    np.savetxt(O/'band_grid_checks.txt',ints,fmt='%.17e',header='Diagnostic band integrals only; no production interpolation.\nColumns: lambda_lo_nm lambda_hi_nm SRW_baseline_J_per_sr SRW_dense_J_per_sr MATLAB_LW_baseline_J_per_sr MATLAB_LW_dense_J_per_sr SRW_baseline_div_dense_minus1 MATLAB_baseline_div_dense_minus1')
    assert np.array_equal(np.loadtxt(O/'band_grid_checks.txt'),np.asarray(ints))
    # Refine points where either pointwise benchmark residual or integration discrepancy is largest.
    ids=np.unique(np.r_[np.arange(0,20,2),np.linspace(0,4999,24).round().astype(int),
        np.argsort(abs(W/Wc-1))[-20:],np.argsort(abs(W/new-1))[-20:]])
    np.save(O/'refinement_indices.npy',ids);np.save(O/'lambda_refinement.npy',lam[ids])
    summary['refinement_nodes']=len(ids)
    (O/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k!='bands'},indent=2),flush=True)

if __name__=='__main__':main()
