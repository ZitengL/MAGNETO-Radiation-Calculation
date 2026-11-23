"""Preview of the complete original grid while independent controls continue."""
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parent/'.matplotlib'))
import analysis as a
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
l,f,_=a.collect('baseline',.234375e-6);w=(abs(f)**2*a.m.HBAR*1e9*121).sum(0)
r=np.load(a.O/'reference_original.npz')['energy_xy']
v=loadmat(a.O/'matlab_baseline_ct1.875um.mat',simplify_cells=True)['result']['W_total_field'][:2].sum(0)
fig,ax=plt.subplots(2,1,figsize=(10,8),sharex=True,layout='constrained')
ax[0].loglog(l,r,label='Original MATLAB, acceleration x+y',lw=1.5)
ax[0].loglog(l,v,'k--',label='Added MATLAB, acceleration+velocity x+y',lw=1)
ax[0].loglog(l,w,label='Native SRW x+y',lw=.8,color='#d95f02')
ax[0].set_ylabel('Energy density [J/(rad/s)/sr]');ax[0].legend(fontsize=8)
ax[1].semilogx(l,(w/r-1)*100,label='SRW / original MATLAB - 1',lw=.65)
ax[1].semilogx(l,(w/v-1)*100,'k',label='SRW / acceleration+velocity - 1',lw=.65)
ax[1].set_ylabel('Pointwise difference [%]');ax[1].set_xlabel('Wavelength [nm]');ax[1].legend(fontsize=8)
for aa in ax:aa.grid(alpha=.25);aa.set_xlim(15,150)
fig.suptitle('PREVIEW: 240 MeV dipole, original5000 wavelengths complete\nIndependent numerical checks still in progress; no fitted correction')
fig.savefig(a.O/'PREVIEW_baseline_comparison.png',dpi=170)
print('Preview written.')
