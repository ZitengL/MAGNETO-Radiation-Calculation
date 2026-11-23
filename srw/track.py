"""Stage 1 only: independent SRW tracking of the two AAC reference cases.

No radiation calculation. Original files on D: are read-only inputs.
Run with the workspace Python 3.12 runtime; see README.md.
"""
from __future__ import annotations
import argparse
from array import array
import hashlib
import json
from pathlib import Path
import struct
import sys
import time
import os

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
SOURCE = Path(__file__).resolve().parents[1]
SRW = ROOT / 'SRW-master/SRW-master'
sys.path.insert(0, str(ROOT / '.srw_dipole_deps'))
sys.path.insert(0, str(SRW / 'env/python'))
os.environ.setdefault('MPLCONFIGDIR', str(OUT / '.matplotlib'))
import numpy as np
from srwpy.srwlib import SRWLParticle, SRWLPrtTrj, SRWLMagFld3D, SRWLMagFldC, srwl

C = 299792458.0
CASES = {'undulator': 'E7GeV_und15p_K1', 'dipole': 'E7GeV_dipole_z1m'}
COMPONENTS = {'x_m': 'arX', 'y_m': 'arY', 'z_m': 'arZ',
              'beta_x': 'arXp', 'beta_y': 'arYp', 'beta_z': 'arZp',
              'Bx_T': 'arBx', 'By_T': 'arBy', 'Bz_T': 'arBz'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def gdf_blocks(path):
    with Path(path).open('rb') as f:
        head = f.read(48)
        assert len(head) == 48 and struct.unpack('<I', head[:4])[0] == 94325877
        while raw := f.read(24):
            assert len(raw) == 24
            name, kind, size = struct.unpack('<16sII', raw)
            payload = f.read(size)
            assert len(payload) == size
            yield name.split(b'\0')[0].decode('ascii'), kind, payload


def reference(case):
    folder = SOURCE / 'gpt_sims/single_electron' / CASES[case]
    stats = {}
    for name, kind, raw in gdf_blocks(folder / 'stats.gdf'):
        if kind & 0x800 and kind & 255 == 3:
            assert name not in stats
            stats[name] = np.frombuffer(raw, dtype='<f8').copy()
    assert np.all(np.diff(stats['time']) > 0)
    p = {}
    for name, kind, raw in gdf_blocks(folder / 'result.gdf'):
        if name in ('m', 'q', 'nmacro') and kind & 0x800 and kind & 255 == 3:
            p[name] = float(np.frombuffer(raw, dtype='<f8')[0])
        if len(p) == 3:
            break
    assert p['nmacro'] == 1 and p['q'] < 0
    p['gamma'] = float(stats['avgG'][0])
    p['gamma_range'] = float(np.ptp(stats['avgG']))
    p['ct_end_m'] = float(C * stats['time'][-1])
    p['baseline_points'] = len(stats['time'])
    return stats, p


def prepare():
    OUT.mkdir(exist_ok=True)
    src = SOURCE / '3d_field_map_of_15_periods_undulator.txt'
    data = np.loadtxt(src, skiprows=1)
    assert data.shape == (31 * 21 * 3001, 7)
    axes = [np.unique(data[:, j]) * 1e-3 for j in range(3)]
    nx, ny, nz = map(len, axes)
    coords = data[:, :3].reshape(nx, ny, nz, 3) * 1e-3
    for j, ax in enumerate(axes):
        shape = [1, 1, 1]; shape[j] = len(ax)
        assert np.allclose(coords[..., j], ax.reshape(shape), rtol=0, atol=1e-15)
    b = data[:, 3:6].reshape(nx, ny, nz, 3)
    # Exact active rotation +90 degrees about z, matching GPT's ECS:
    # (xw,yw,zw)=(-yl,xl,zl+1); (Bxw,Byw,Bzw)=(-Byl,Bxl,Bzl).
    w = b.transpose(2, 0, 1, 3)[:, :, ::-1, :].copy()
    w = w[..., [1, 0, 2]]
    w[..., 0] *= -1
    x, y, z = -axes[1][::-1], axes[0], axes[2]
    np.savez_compressed(OUT / 'undulator_field_rotated.npz', x=x, y=y, z=z, B=w)
    refs = {}
    for case in CASES:
        s, p = reference(case)
        np.savez_compressed(OUT / f'{case}_gpt_reference.npz', **s)
        refs[case] = p
    provenance = {'inputs': {str(p): sha(p) for p in [src,
        SOURCE / 'field_map/und15p_K1_field_map.gdf',
        Path(os.environ.get('GPT_ROOT', 'GPT'), r'elems\sources\BStatic\rectmagnet.c'),
        Path(os.environ.get('GPT_ROOT', 'GPT'), r'kernel\gps.h'),
        SRW / 'env/python/srwpy/srwlpy.pyd',
        SRW / 'cpp/src/core/srgtrjdt.cpp']}, 'reference': refs,
        'rotation': 'xw=-yl; yw=xl; zw=zl+1; Bxw=-Byl; Byw=Bxl; Bzw=Bzl',
        'field_shape_zyx': list(w.shape),
        'axes_max_uniform_grid_deviation_m': [float(np.max(np.abs(a-np.linspace(a[0],a[-1],len(a))))) for a in (x,y,z)],
        'constants_policy': 'Exact GPT gamma, native SRW electron charge/mass, relE0=1; no fitted field or energy correction.'}
    (OUT / 'inputs.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    print(json.dumps(provenance, indent=2), flush=True)


def field_container(x, y, z, b):
    f = SRWLMagFld3D(*[array('d', np.ascontiguousarray(b[..., j]).ravel()) for j in range(3)],
        _nx=len(x), _ny=len(y), _nz=len(z), _rx=float(np.ptp(x)),
        _ry=float(np.ptp(y)), _rz=float(np.ptp(z)), _nRep=1, _interp=1,
        _arX=array('d', x), _arY=array('d', y), _arZ=array('d', z))
    return SRWLMagFldC([f], array('d', [0]), array('d', [0]), array('d', [1.0]))


def dipole_field(zstep_um, params):
    # Original formula for this case's non-overlapping GPT rectmagnet.
    # Local (X,Y,Z)=(global y,-global x,global z-1), b2=dl=0.
    # Table includes x=0 exactly; the exact symmetric orbit remains in x=0.
    # y grid lies inside the unchanged GPT width |y|<=0.1 m; B is y-independent there.
    x = np.array([-.001, 0., .001])
    y = np.array([-.1, 0., .1])
    z = np.linspace(-.75, .75, round(1.5 / (zstep_um * 1e-6)) + 1)
    a = np.arcsin(.0254 / 11)
    theta = np.arctan(11 * np.sin(a) / (11 * np.cos(a) - 1))
    gamma = params['gamma']
    B0 = np.sin(theta) * np.sqrt(gamma*gamma-1) * params['m'] * C / abs(params['q'])
    d = np.abs(z) - .5
    ef = np.exp(40 * d)[:, None, None]
    h = -40*x[None, None, :]
    bn = B0 / (1+2*ef*np.cos(h)+ef*ef)
    b = np.zeros((len(z), len(y), len(x), 3))
    b[..., 0] = -bn*(1+ef*np.cos(h))
    b[..., 2] = -bn*ef*np.sin(h)*np.where(z[:, None, None]<0, -1., 1.)
    return field_container(x, y, z, b), {'B0_T': float(B0), 'target_angle_rad':float(theta),
        'zstep_um': float(zstep_um), 'table_shape_zyx':list(b.shape),
        'formula': 'GPT rectmagnet_sim; dl=b2=0, b1=40/m; exact cutoff at |z-1|>0.75 m',
        'transverse_domain': 'x=+-1 mm contains the exact symmetry plane; y=+-0.1 m equals GPT width; no transverse averaging'}


def run(case, factor, fieldstep):
    provenance = json.loads((OUT/'inputs.json').read_text())
    p = provenance['reference'][case]
    label = f'{case}_ct{60/factor:g}um' + (f'_field{fieldstep:g}um' if case=='dipole' else '')
    t0=time.perf_counter()
    if case == 'undulator':
        f = np.load(OUT/'undulator_field_rotated.npz')
        mag=field_container(f['x'],f['y'],f['z'],f['B'])
        info={'source': 'complete original 3D field rotated +90 degrees; trilinear interpolation'}
    else:
        mag,info=dipole_field(fieldstep,p)
    tr=SRWLPrtTrj(_partInitCond=SRWLParticle(_gamma=p['gamma'],_relE0=1,_nq=-1))
    n=(p['baseline_points']-1)*factor+1
    tr.allocate(n,True)
    tr.ctStart=0.; tr.ctEnd=p['ct_end_m']
    srwl.CalcPartTraj(tr,mag,[1]) # SRW manual fixed-step RK4, no user integrator
    a={'ct_m':np.linspace(tr.ctStart,tr.ctEnd,n)}
    a.update({k:np.asarray(getattr(tr,v)).copy() for k,v in COMPONENTS.items()})
    assert all(np.all(np.isfinite(v)) for v in a.values())
    np.savez_compressed(OUT/f'{label}.npz',**a)
    metrics={'case':case,'factor':factor,'points':n,'ct_step_um':float((tr.ctEnd-tr.ctStart)/(n-1)*1e6),
        'seconds':time.perf_counter()-t0,'gamma':p['gamma'],'field':info,
        'final':{k:float(v[-1]) for k,v in a.items()},
        'ranges':{k:[float(v.min()),float(v.max())] for k,v in a.items()},
        'max_abs_beta_norm_error':float(np.max(np.abs(a['beta_x']**2+a['beta_y']**2+a['beta_z']**2-(1-p['gamma']**-2))))}
    (OUT/f'{label}.json').write_text(json.dumps(metrics,indent=2))
    print(json.dumps(metrics,indent=2),flush=True)


def differences(a,b,thin=1):
    return {k:float(np.max(np.abs(a[k]-b[k][::thin]))) for k in COMPONENTS if not k.endswith('_T')}


def report():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    summary={}
    for case in CASES:
        files=sorted(OUT.glob(f'{case}_ct*.json'))
        if not files: continue
        info=[json.loads(f.read_text()) for f in files]
        selected=[i for i in info if case=='undulator' or i['field']['zstep_um']==15]
        selected.sort(key=lambda v:v['factor'])
        labels=[f"{case}_ct{60/i['factor']:g}um"+(f"_field{i['field']['zstep_um']:g}um" if case=='dipole' else '') for i in selected]
        arrays=[dict(np.load(OUT/f'{l}.npz')) for l in labels]
        finest=arrays[-1]; fineinfo=selected[-1]
        g=dict(np.load(OUT/f'{case}_gpt_reference.npz'))
        ref={'x_m':g['avgx'],'y_m':g['avgy'],'z_m':g['avgz'],
             'beta_x':g['avgBx'],'beta_y':g['avgBy'],'beta_z':g['avgBz']}
        gf={k:np.interp(finest['ct_m'],C*g['time'],v) for k,v in ref.items()}
        rows=[]
        for i,a in zip(selected,arrays):
            rows.append({'ct_step_um':i['ct_step_um'],'seconds':i['seconds'],
                'vs_finest_max_absolute':differences(a,finest,fineinfo['factor']//i['factor'])})
        shared={k:finest[k][::fineinfo['factor']] for k in ref}
        dg={k:float(np.max(np.abs(shared[k]-ref[k]))) for k in ref}
        row={'runs':info,'tracking_convergence_fixed_field':rows,'finest_vs_GPT_at_original_times':dg}
        if case=='dipole':
            fieldruns=sorted([i for i in info if i['factor']==4],key=lambda i:-i['field']['zstep_um'])
            row['field_convergence_fixed_ct15um']=[]
            for i in fieldruns:
                a=dict(np.load(OUT/f"dipole_ct15um_field{i['field']['zstep_um']:g}um.npz"))
                row['field_convergence_fixed_ct15um'].append({'field_step_um':i['field']['zstep_um'],
                    'vs_field15um_max_absolute':differences(a,finest)})
        summary[case]=row
        coord='x' if case=='undulator' else 'y'
        field='By_T' if case=='undulator' else 'Bx_T'
        scale=1e6 if case=='undulator' else 1e3
        unit='um' if case=='undulator' else 'mm'
        fig,axes=plt.subplots(3,2,figsize=(12,10),layout='constrained')
        z=finest['z_m']
        panels=[(field,1.,'Magnetic field [T]'),(coord+'_m',scale,f'Orbit {coord} [{unit}]'),('beta_'+coord,1e3,r'$\beta_'+coord+r'$ [mrad]')]
        for j,(key,s,ylabel) in enumerate(panels):
            gkey={'Bx_T':'avgfBx','By_T':'avgfBy'}.get(key)
            gv=g[gkey] if gkey else ref[key]
            axes[j,0].plot(g['avgz'],gv*s,color='0.5',lw=2,label='GPT reference')
            axes[j,0].plot(z,finest[key]*s,'--',lw=1,label='SRW, ct step 15 um')
            axes[j,0].set_ylabel(ylabel)
            if j==0:
                d=finest[key][::fineinfo['factor']]-gv
                axes[j,1].plot(g['avgz'],d,lw=1)
                axes[j,1].set_ylabel('SRW - GPT field [T]')
            else:
                for i,a in zip(selected[:-1],arrays[:-1]):
                    thin=fineinfo['factor']//i['factor']
                    axes[j,1].plot(a['z_m'],(a[key]-finest[key][::thin])*(1e9 if j==1 else 1e6),label=f"ct {60/i['factor']:g} - 15 um")
                axes[j,1].set_ylabel('Orbit difference [nm]' if j==1 else r'$\Delta\beta$ [$\mu$rad]')
        for ax in axes.flat:
            ax.grid(alpha=.25); ax.set_xlabel('Global z [m]')
        axes[0,0].legend(); axes[1,1].legend()
        fig.suptitle(f'{case.capitalize()}: independent SRW particle tracking (no radiation)')
        fig.savefig(OUT/f'{case}_trajectory.png',dpi=170); plt.close(fig)
        np.savetxt(OUT/f'{case}_trajectory.csv',np.column_stack(list(finest.values())),delimiter=',',header=','.join(finest),comments='',fmt='%.16e')
    (OUT/'tracking_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:{'tracking_convergence':v['tracking_convergence_fixed_field'],'vs_GPT':v['finest_vs_GPT_at_original_times']} for k,v in summary.items()},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('action',choices=['prepare','run','report'])
    parser.add_argument('--case',choices=list(CASES),default='undulator')
    parser.add_argument('--factor',type=int,default=1)
    parser.add_argument('--fieldstep',type=float,default=15)
    args=parser.parse_args()
    if args.action=='prepare': prepare()
    elif args.action=='run': run(args.case,args.factor,args.fieldstep)
    else: report()
