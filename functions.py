import matplotlib
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit
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
import time
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from scipy.stats import norm
import matplotlib.mlab as mlab
import sys
from astropy.table import hstack

#--------------------------------------------------------------------------
#Functions

def multivar_gauss_error(nnsol,A):
    covmask=nnsol==0
    covmask2=nnsol>0
    nnsol_mask=nnsol[covmask2]
    C=A[:,covmask2]
    Mdust_cov=[]
    Mgas_cov=[]
    Lir_cov=[]
    chi2_cov=[]
    eLir=-99
    eMD=-99
    eMG=-99
    max_allowed = 500
    attempt = 0
    loop_time_start=time.time()
    while (eLir<1.) or np.isnan(eLir) or (eMD<1.) or np.isnan(eMD):
        attempt+=1
        if attempt>max_allowed or time.time()-loop_time_start>1200:
            print('Max attempts exceeded: TIMEOUT ERROR')
            break


        try:
            cov = np.matrix(np.dot(C.T, C)).I.A
            covsam=np.random.multivariate_normal(nnsol_mask,cov,attempt*10**3)
            incl_val=[]
            for val,_ in enumerate(covsam[:,0]):
                if all(h>0 for h in covsam[val,:]):
                    incl_val.append(val)
            covsam=covsam[incl_val,:]
            nnsol_cov=np.zeros((len(covsam[:,0]),len(nnsol)))
            nnsol_cov[:,covmask2]=covsam
        except:
            nnsol_cov=nnsol
            break

        try:
            nnsol_cov_sum=np.sum(nnsol_cov[:,total-irtemp:],axis=1)
            Mdust_cov=DGratio*nnsol_cov_sum*b*(const.m_p/const.M_sun)
            for x,obj in enumerate(nnsol_cov[:,0]):
                Lir_cov.append(LumIR_(RFull,full_template(*nnsol_cov[x,:]),z,cosmo))
                chi2_cov.append(np.sum((f_opt(wav_ar,*nnsol_cov[x,:])-flux)**2*(flux_e**2)**-1))
            ConfIntL=mean_confidence_interval_fast(Lir_cov/Lir_nnls)
            ConfIntM=mean_confidence_interval_fast(Mdust_cov/Mdust_nnls)
            eLir=np.mean([ConfIntL[0]-ConfIntL[1],ConfIntL[2]-ConfIntL[0]])*np.median(Lir_cov)
            eMD=np.mean([ConfIntM[0]-ConfIntM[1],ConfIntM[2]-ConfIntM[0]])*np.median(Mdust_cov)
            eMG=eMD*deltaGDR

        except:
            try:
                gaussM=norm.fit(Mdust_cov/Mdust_nnls)
                gaussL=norm.fit(Lir_cov/Lir_nnls)
                eLir=gaussL[1]*np.median(Lir_cov/Lir_nnls)
                eMD=gaussM[1]*np.median(Mdust_cov/Mdust_nnls)
                eMG=deltaGDR*eMD
            except:
                print('FAILED')
                break
    return eLir,eMD,eMG,nnsol_cov,attempt




def convolver(sedx,sedy,lam,T,obs):
    sedx=np.array(sedx)
    sedy=np.array(sedy)
    lam=np.array(lam)
    T=np.array(T)
    lam=lam[T>0.0]
    T=T[T>0.0]
    if len(sedx)!=len(lam):
        sedy=np.interp(lam,sedx,sedy)

        #sedy=f(lam)
    #integrator=sp.integrate.simps
    integrator=np.trapz
    #return obs*integrator(T*sedy*lam**-1,lam)/integrator(T,lam)
    return integrator(T*sedy/lam,lam)/integrator(T/lam,lam)


def integrate_filter(sedx,sedy,lam,T,scale=1., z=0):
        """
        Integrate the template through a `FilterDefinition` filter object.

        The `grizli` interpolation module should be used if possible:
        https://github.com/gbrammer/grizli/
        """
        try:
            import grizli.utils_c
            interp = grizli.utils_c.interp.interp_conserve_c
        except ImportError:
            interp = grizli.utils.interp_conserve

        templ_filter = interp(lam, sedx*(1+z),
                              sedy*scale, left=0, right=0)

        # f_nu/lam dlam == f_nu d (ln nu)
        integrator = np.trapz
        temp_int = integrator(T*templ_filter/lam, lam)/integrator(T/lam, lam)

        return temp_int


def interp_conserve(x, xp, fp, left=0., right=0.):
    midpoint = (x[1:]-x[:-1])/2.+x[:-1]
    midpoint = np.append(midpoint, np.array([x[0],x[-1]]))
    midpoint = midpoint[np.argsort(midpoint)]
    int_midpoint = np.interp(midpoint, xp, fp, left=left, right=right)
    int_midpoint[midpoint > xp.max()] = right
    int_midpoint[midpoint < xp.min()] = left

    fullx = np.append(xp, midpoint)
    fully = np.append(fp, int_midpoint)

    so = np.argsort(fullx)
    fullx, fully = fullx[so], fully[so]

    outy = x*0.
    dx = midpoint[1:]-midpoint[:-1]
    for i in range(len(x)):
        bin = (fullx >= midpoint[i]) & (fullx <= midpoint[i+1])
        outy[i] = np.trapz(fully[bin], fullx[bin])/dx[i]

    return outy


#import functions_c

def convolver_(sedx,sedy,lam,T,obs):
    interp=interp_conserve
    templ_filter = interp(lam, sedx,sedy, left=0, right=0)
    integrator = np.trapz
    temp_int = integrator(T*templ_filter/lam, lam)/integrator(T/lam, lam)

    return temp_int


def find_nearest(array, value):
    n = [abs(i-value) for i in array]
    idx = n.index(min(n))
    return idx

def mean_confidence_interval(data, confidence=0.68):
    a = np.sort(1.0 * np.array(data))
    n = len(a)
    m, se = np.median(a), sp.stats.sem(a)
    dlo = a[a < m]
    dhi=a[a>m]
    frac = 0
    count = 0
    while frac < 1-confidence:
        frac = len(dlo[:count]) / len(dlo)
        count += 1
        hlo=a[np.argwhere(a==dlo[count])]
    hlo=float(hlo[0])
    frac = 0
    count = 0
    while frac < confidence:
        frac = len(dhi[:count]) / len(dhi)
        count += 1
        hhi=a[np.argwhere(a==dhi[count])]
    hhi=float(hhi[0])
    return m, hlo, hhi

def mean_confidence_interval_fast(data, confidence=0.68):
    a = np.sort(1.0 * np.array(data))
    n = len(a)
    m=np.median(a)
    lo=0.5-confidence/2
    hi=0.5+confidence/2
    hlo=a[int(lo*len(a))]
    hhi=a[int(hi*len(a))]
    return m, hlo, hhi

def mean_confidence_interval_fast_error(sols,zero_sol, confidence=0.68,mode='new'):
    if mode=='legacy':
        data=np.array(sols)/zero_sol
    if mode=='new':
        data=np.array(sols)/np.median(sols)
    a = np.sort(1.0 * np.array(data))
    n = len(a)
    m=np.median(a)
    lo=0.5-confidence/2
    hi=0.5+confidence/2
    hlo=a[int(lo*len(a))]
    hhi=a[int(hi*len(a))]
    err=np.mean([m-hlo,hhi-m])*np.median(sols)
    return err

def mgas_metallicity_manucci(Mstar,Mdust,SFR):
    if Mstar>0:
        mu_032=np.log10(Mstar)-0.32*np.log10(SFR)
        x=mu_032-10.
        if mu_032<10.2:
            met_KD02=8.90+0.47*x
        if mu_032>10.5:
            met_KD02=9.07
        else:
            met_KD02=8.90+0.39*x-0.2*x**2-0.077*x**3+0.064*x**4

        #Convert from KD02 to PP04 N2 using eq.1 and table 3 from Kewley+08
        met_PP04=569.4927-192.51820*met_KD02+21.918360*met_KD02**2-0.8278840*met_KD02**3
        deltaGDR=10**(10.54-0.99*met_PP04)
        try:
            Mgas=Mdust*deltaGDR
        except:
            Mgas=-99
    else:
        Mgas=-99
        deltaGDR=-99
    return Mgas,deltaGDR

def LumIR_(sedx,sedy,z,cosmology):
    #Computes TIR of an SED
    #sedx - rest frame wavelength
    #sedy - sed values in mJy
    #z - redshift of the source
    #cosmo - cosmology
    flux=sedy*10**-26
    lum_dist=(cosmology.luminosity_distance(z)).to(u.cm).value
    luminosity_v=((4*np.pi*lum_dist**2*flux*(1+z)**-1)*u.erg).si
    wavrange=[find_nearest(sedx,8),find_nearest(sedx,1000)+1]
    luminosity_v=luminosity_v[wavrange[0]:wavrange[1]]
    wavelength=((sedx[wavrange[0]:wavrange[1]]*(1.+0))*u.um).si
    Luminosity_IR=integrate.trapz((luminosity_v*const.c*wavelength**-2).value,wavelength.value)*(3.839*10**26)**-1

    return Luminosity_IR

def SFRMS_Schreiber_S(z,Ms): #Schreiber Main Sequence input Chabrier output Chabrier
    Ms=0.61**-1*Ms
    r=np.log10(1+z)
    m=np.log10(Ms/10**9)
    m0=0.5
    a0=1.5
    a1=0.3
    m1=0.36
    a2=2.5
    q=m-m1-a2*r
    if q<0:
        q=0
    logSFR=(m-m0)+a0*r-a1*(q)**2
    return 0.63*(10**logSFR) #Output sfr

def SFRMS_Schreiber(z,Ms): #Schreiber Main Sequence input Chabrier output Chabrier
    Ms=0.61**-1*np.array(Ms)
    r=np.log10(1+z)
    m=np.log10(Ms/10**9)
    m0=0.5
    a0=1.5
    a1=0.3
    m1=0.36
    a2=2.5
    q=m-m1-a2*r
    q[q<0]=0
    logSFR=(m-m0)+a0*r-a1*(q)**2
    return 0.63*(10**logSFR) #Output sfr


def radio_slope_delhz(z,LIR,dl,alpha):
    dl=(dl*u.cm).to(u.m)
    q=2.88*(1+z)**(-0.19)
    L14=(10**(np.log10(LIR*(3.839*10**26)/(3.75*10**12))-q))*u.W/u.Hz
    S3=L14*(1+z)**(alpha+1)*(4*np.pi*dl**2)**-1*(3/1.4)**alpha
    S3=S3.to(u.mJy).value
    A=S3*(10**5)**alpha
    return A

def radio_slope_delv(z,LIR,Mstar,dl,alpha):
    dl=(dl*u.cm).to(u.m)
    A=(1+z)
    B=np.log10(Mstar)-10
    q=2.646*A**(-0.023)-B*0.148
    L14=(10**(np.log10(LIR*(3.839*10**26)/(3.75*10**12))-q))*u.W/u.Hz
    S3=L14*(1+z)**(alpha+1)*(4*np.pi*dl**2)**-1*(3/1.4)**alpha
    S3=S3.to(u.mJy).value
    A=S3*(10**5)**alpha
    return A
'''
def rad_flux(lam,method='delv20',alpha=-0.75,use_K=False):

    #Calculate the radio contribution:
    #Two methods available:
    #Delvecchio+20, with Lir and Mstar dependence - delv20
    #Delhaize+17, only with Lir dependence - delhz17
    #If no stellar mass available, can use the K-band prediction instead - use_K

    if method=='delv20':
        slope=radio_slope_delv
        if Mstar<0:
            if use_K:
                output=slope(z,Lir_total,mass_K,lum_dist,alpha)*lam**(-alpha)
            else:
                print('Cannot use Delvecchio+20 radio method, without M*')
                output=lam*0.
        else:
            output=slope(z,Lir_total,Mstar,lum_dist,alpha)*lam**(-alpha)

    elif method=='delhz17':
        slope=radio_slope_delhz
        output=slope(z,Lir_total,lum_dist,alpha)*lam**(-alpha)
    return output
'''

def chi_vectors(A,solutions,N=10**3):
    vectors=[]
    for i,solution in enumerate(solutions[0,:]):
        solution=solutions[:,i]
        covmask=np.argwhere(solution==0)
        covmask2=solution>0
        nnsol_mask=np.delete(solution,covmask,None)
        C=np.delete(A,covmask,axis=1)
        try:
            cov = np.matrix(np.dot(C.T, C)).I.A
            covsam=np.random.multivariate_normal(nnsol_mask,cov,N)
            incl_val=[]
            for val in range(len(covsam[:,0])):
                if all(h>0 for h in covsam[val,:]):
                    incl_val.append(val)
            covsam=covsam[incl_val,:]
            nnsol_cov=np.zeros((len(covsam[:,0]),len(solution)))
            nnsol_cov[:,covmask2]=covsam

        except Exception as e:
            print(e)
            nnsol_cov=solution

        vectors.append(nnsol_cov)
    stack=np.vstack(vectors)
    return stack


def f_pdr(gamma,Umin):
    return (gamma*np.log(10**6/10**2))*((1-gamma)*(1-Umin/10**6)+gamma*np.log(10**6/Umin))**-1

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def moving_average2(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def filter_cont(lmbd,w,dx):
    lmbd=np.array(lmbd)
    fltrx=[]
    fltry=[]
    for i,_ in enumerate(lmbd):
        c=lmbd[i]
        fltr = np.arange(c-w, c+w, dx)
        fltrtop = np.ones_like(fltr)
        fltrtop[(fltr<c-w) | (fltr>c+w)]=0
        fltrx.append(fltr)
        fltry.append(fltrtop)

    fltrx=np.array(fltrx)
    fltry=np.array(fltry)
    of=np.array([fltrx,fltry])
    of=np.swapaxes(of,0,1)
    return Table(of)



def make_output_table(table_out):
    outputnames=['ID','LIR_total','eLIR_total','MD',
                 'eMD','z','chi2_red','f_agn','efagn','lastdet','MG','eMG',
                 'deltaGDR','attempts','Mstar','fgas','fgas_FMR','Lir_med',
                 'eLir68','Mdust_med','eMdust68','Umin','qpah','gamma','U','sU',
                 'Lagn','eLagn','Lir_draine','eLir_draine','mass','e_mass','mass_K','chi2','nfilt','sfr_opt','Av']

    dtype=[float]*len(outputnames)
    dtype[0]=int

    table_0=Table(np.zeros(len(outputnames)), names=outputnames,dtype=dtype) #Make empty table
    table_0.write(table_out,overwrite=True)
    return None

def get_qso_templates():
    qso_templ1=Table.read('templates/QSO/shen2016_ext_Av0.0.fits')
    qso_templ2=Table.read('templates/QSO/shen2016_ext_Av0.1.fits')
    qso_templ3=Table.read('templates/QSO/shen2016_ext_Av0.2.fits')
    qso_templ4=Table.read('templates/QSO/shen2016_ext_Av0.3.fits')
    qso_templ5=Table.read('templates/QSO/shen2016_ext_Av0.4.fits')

    qso_wave=qso_templ1['wave']

    QSO=np.zeros((5,len(qso_wave)))

    QSO[0,:]=qso_templ1['flux']
    QSO[1,:]=qso_templ2['flux']
    QSO[2,:]=qso_templ2['flux']
    QSO[3,:]=qso_templ3['flux']
    QSO[0,:]=qso_templ4['flux']
    return QSO,qso_wave


def convert_to_mjy(unit):
    if unit=='mJy':
        conv=1.
    elif unit=='uJy':
        conv=(u.uJy).to(u.mJy)
    elif unit=='nJy':
        conv=(u.nJy).to(u.mJy)
    else:
        sys.exit(f'Unrecognised unit - {unit}, accepted formats: uJy, mJy, nJy')
    return conv
