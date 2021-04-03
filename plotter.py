from astropy.io import fits
from astropy.table import Table
import numpy as np
import scipy as sc
from scipy import stats
from scipy.stats import norm
from astropy  import units as u
from astropy import constants as const
import scipy as sp
from scipy.integrate import quad
from scipy import integrate
import math
from scipy.signal import convolve
from scipy import interpolate
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from scipy.stats import norm
from matplotlib.pyplot import *
import  matplotlib.pyplot as plt

from init import *
from config import *

from filter_importer import *


for i,_ in enumerate(P[:,0]):
    try:
        sed_file=Table.read(f'{sedloc}{int(_)}.fits',memmap=True)
    except:
        print(f'ID{_} not found')
        continue

    fig= plt.figure(figsize=(10,5))
    ax = fig.add_subplot(1,1,1)
    textypos=0.9
    textxpos=0.6
    textsep=0.08


    wave=sed_file['lambda']

    SF=sed_file['stellar']
    AGN=sed_file['AGN']
    IR=sed_file['IR']
    TOTAL=sed_file['Total']+sed_file['Radio']

    ax.fill_between(wave,0,SF,color='royalblue',alpha=0.2,label='Stellar')
    ax.fill_between(wave,0,AGN,color='g',alpha=0.2,label='AGN')
    ax.fill_between(wave,0,IR,color='maroon',alpha=0.2,label='Dust')
    ax.plot(wave,TOTAL,'k',label='Total',lw=3,alpha=0.6,zorder=10)

    ax.text(0.02, 0.85,'ID '+str(int(_)), color='k',fontsize=20,transform=ax.transAxes)
    ax.text(0.02, 0.75,r'$z$={:.2f}'.format(P[i,1]), color='k',fontsize=20,transform=ax.transAxes)

    flux=G[i,:]
    flux_e=E[i,:]
    sfx=W[i,:]
    points=((flux/flux_e)>=3)


    ax.errorbar(sfx[points],flux[points],yerr=flux_e[points],color='red',fmt='s',capsize=5,capthick=1,ms=12,markerfacecolor='white',mew=2,barsabove=True)
    ax.scatter(sfx[~points],(flux+3*flux_e)[~points],marker=r'$\downarrow$',s=300,color='red',zorder=11)


    ax.set_ylabel(r'$f_{\nu}$ [mJy]',fontsize=25)
    ax.set_xlabel(r'$\lambda_{obs}$ $[\mu m]$',fontsize=25)
    ax.legend(fontsize=12)
    #ax.set_ylim(10**-4,10**3)
    #ax.set_xlim(.5,10**5.7)
    ax.set_ylim(10**-5,10**3)
    ax.set_xlim(.1,10**5.7)
    #ax.set_xlim(0.9,1e4)
    #ax.set_ylim(10**-3,10)
    ax.grid(alpha=0.4)
    ax.set_yscale('log')
    ax.set_xscale('log')

    plt.tight_layout()
    plt.savefig(f'{figloc}/{_}.pdf')
    plt.clf()
