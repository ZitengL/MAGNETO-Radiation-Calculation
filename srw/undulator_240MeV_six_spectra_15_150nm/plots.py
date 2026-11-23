"""Static paper-oriented plots of direct native SRW spectra and raw comparisons."""
import calculate as m
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import maximum_filter1d
O=m.O

def save(fig,name):
    fig.savefig(O/(name+'.png'),dpi=190);fig.savefig(O/(name+'.svg'));plt.close(fig)

def draw(data,ref,stats):
    lam=data['lambda_nm'];S=data['SRW_energy_xy'];M=data['MATLAB_xy'];T=data['MATLAB_total'];ctrl=data['SRW_integration_control']
    colors=plt.get_cmap('tab10').colors
    fig,ax=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for i in range(6):
        for a in ax:a.plot(lam,S[i],color=colors[i],lw=1.1,label=f'{i/10:.1f} mrad')
    for a in ax:a.set_xscale('log');a.set_xlim(15,150);a.set_xlabel('Wavelength [nm]');a.grid(alpha=.25);a.set_ylabel(r'$d^2W/(d\omega\,d\Omega)$ [J/(rad/s)/sr]')
    ax[0].set_yscale('log');ax[0].legend(fontsize=8);ax[1].legend(fontsize=8)
    fig.suptitle('240 MeV, 15-period undulator: six native SRW spectra\nCompensated gamma and magnetic-field inputs, observers on R = 11 m sphere')
    save(fig,'spectra_overlay')
    fig,axes=plt.subplots(3,2,figsize=(13,11),layout='constrained',sharex=True)
    for i,a in enumerate(axes.flat):
        a.loglog(lam,S[i],lw=1.4,label='Native SRW x+y')
        a.loglog(lam,M[i],'--',lw=1,label='Original MATLAB x+y')
        a.loglog(lam,T[i],':',lw=.8,label='Original MATLAB total x+y+z')
        a.set_title(f'theta_x = {i/10:.1f} mrad');a.grid(alpha=.25);a.set_xlim(15,150)
        a.set_ylabel('Energy density [J/(rad/s)/sr]');a.set_xlabel('Wavelength [nm]')
        if i==0:a.legend(fontsize=7)
    fig.suptitle('240 MeV undulator: direct SRW / MATLAB comparison\nSame wavelength nodes and exact spherical observer geometry; no fitted shift or scaling')
    save(fig,'six_angles_MATLAB_comparison')
    fig,axes=plt.subplots(3,2,figsize=(13,10),layout='constrained',sharex=True)
    half=round(.01/np.mean(np.diff(np.log(lam))))
    for i,a in enumerate(axes.flat):
        a.semilogx(lam,1e6*(S[i]-M[i])/max(M[i]),lw=.85,label='SRW - MATLAB x+y')
        a.semilogx(lam,1e6*(S[i]-T[i])/max(M[i]),'--',lw=.65,label='SRW - MATLAB total')
        a.set_title(f'theta_x = {i/10:.1f} mrad');a.set_ylabel('Difference / MATLAB peak [ppm]');a.set_xlabel('Wavelength [nm]');a.grid(alpha=.25);a.set_xlim(15,150)
        if i==0:a.legend(fontsize=7)
    fig.suptitle('Raw spectral differences relative to each MATLAB peak\nNo amplitude or wavelength fit; all original samples retained')
    save(fig,'six_angles_differences')
    fig,axes=plt.subplots(3,2,figsize=(13,10),layout='constrained',sharex=True)
    for i,a in enumerate(axes.flat):
        a.semilogx(lam,100*(S[i]/M[i]-1),lw=.75)
        a.set_title(f'theta_x = {i/10:.1f} mrad');a.set_ylabel('Raw (SRW / MATLAB x+y - 1) [%]')
        a.set_xlabel('Wavelength [nm]');a.grid(alpha=.25);a.set_xlim(15,150)
    fig.suptitle('Raw pointwise relative differences, all original wavelength nodes\nWeak minima have small denominators; no intensity threshold or sample removal')
    save(fig,'pointwise_relative_errors')
    fig,axes=plt.subplots(2,1,figsize=(11,8),layout='constrained',sharex=True)
    for i in range(6):
        env=maximum_filter1d(S[i],size=2*half+1,mode='nearest')
        axes[0].semilogx(lam,1e6*(S[i]-ctrl[i])/env,lw=.7,label=f'{i/10:.1f} mrad')
        me=maximum_filter1d(M[i],size=2*half+1,mode='nearest')
        axes[1].semilogx(lam,100*(S[i]-M[i])/me,lw=.7,label=f'{i/10:.1f} mrad')
    axes[0].set_ylabel('SRW (step30 - step60) / envelope [ppm]')
    axes[1].set_ylabel('(SRW - MATLAB x+y) / envelope [%]');axes[1].set_xlabel('Wavelength [nm]')
    for a in axes:a.grid(alpha=.25);a.set_xlim(15,150);a.legend(fontsize=7,ncol=3)
    fig.suptitle('240 MeV undulator: integration check and benchmark residual\nEnvelope is max within +/-1% in log wavelength; formal spectra are not smoothed')
    save(fig,'convergence_and_residual')
    for i in range(6):
        for log in (True,False):
            fig,a=plt.subplots(figsize=(10,5.5),layout='constrained')
            a.plot(lam,S[i],lw=1.4,label='Native SRW x+y')
            a.plot(lam,M[i],'--',lw=1,label='MATLAB x+y')
            for n in ([1,3,5] if i==0 else [1,2,3,4,5]):
                h=ref['harmonics_nm'][i,n-1]
                if 15<=h<=150:a.axvline(h,color='.65',lw=.6,ls=':')
            a.set_xscale('log')
            if log:a.set_yscale('log')
            a.set_xlim(15,150);a.grid(alpha=.25);a.legend(fontsize=8);a.set_xlabel('Wavelength [nm]')
            a.set_ylabel(r'$d^2W/(d\omega\,d\Omega)$ [J/(rad/s)/sr]')
            a.set_title(f'240 MeV undulator, theta_x={i/10:.1f} mrad, R=11m\nDotted: original MATLAB harmonic markers; no fitted wavelength shift')
            save(fig,f'spectrum_theta_{i/10:.1f}mrad_'+('logy' if log else 'liny'))
