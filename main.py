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
from astropy.table import vstack
from functions import *



timer_start=time.time()

from init import *
from config import *

if save_table:
    make_output_table(table_out)
#--------------------------------------------------------------------------

assert len(sfx)==len(FILTERS)==len(band_names)==len(err_band_names),'ERROR: Filter range does not match available photometry, check band names'


print('Photometry range ok')

#---------------------------------------------------------
#Import templates
from template_importer import *
#---------------------------------------------------------
def galproc(galaxy):
    galaxy_index=galaxy
    if verbose:
        print ('Calculating index:',(P[galaxy_index,0]))



    z=P[galaxy_index,1]
    if z<=0:
        print(f'ID:{P[galaxy_index,0]} has z<=0, skipped')
        return None

    try:
        Mstar=P[galaxy_index,2]
    except:
        Mstar=-99

    if extra_bands:
        sfx_extra=W_ex[galaxy_index,:]

    #Redshift the template range
    RFullred=RFull*(1+z)
    start = time.time()


    #=======================================
    #



    flux_i=np.copy(G[galaxy_index,:])
    flux_i_e=np.copy(E[galaxy_index,:])
    flux_i_e_orig=np.copy(flux_i_e)

    snr=flux_i/flux_i_e
    wav_ar_i=W[galaxy_index,:]

    try:
        lastdet=(wav_ar_i[snr>3][-1])/(1+z)
    except:
        lastdet=-99



    mask=[]
    for i,obj in enumerate(flux_i):
        snr=flux_i[i]/flux_i_e[i]
        if snr<=3. and flux_i[i]>10**-11 and flux_i_e[i]>10**-10:
            flux_i_e[i]=3.*flux_i_e[i] #Upper limit
        elif flux_i[i]<10**-11:
            mask.append(i)
        else:
            continue


    flux=np.delete(flux_i, mask,None)
    flux_e=np.delete(flux_i_e, mask,None)
    flux_e_orig=np.delete(flux_i_e_orig,mask,None)
    wav_ar=np.delete(wav_ar_i,mask,None)

    nfilt=len(flux)


    #return int(P[galaxy_index,0]),nfilt



    #Create extra filters and add them to the existing filter array
    if extra_bands:
        filt_extra=np.array(filter_cont(sfx_extra,5,0.2))
        if verbose==2:
            print(f'Created square waves at {sfx_extra} um')

    #SNR reducer



    for i,obj in enumerate(wav_ar):
        if wav_ar[i]<=25.:
            flux_e[i]=(flux_e[i]**2+(uncert_scale*flux[i])**2)**0.5
            #if flux[i]/flux_e[i]>=5:
            #    flux_e[i]=flux[i]*(1./5.) # REDUCE THE IMPORTANCE OF STELLAR DATA


    lum_dist=(cosmo.luminosity_distance(z)).to(u.cm).value

    try:
        FL=(1.+z)*(4.*pi*lum_dist**2)**-1
    except:
        print('Failed at ID ', P[galaxy_index,0])
        return None



    #==========================================================================
    'Synthetic photometry'
    timer_synphot=time.time()

    igm_factor=RFull*0+1.

    if igm_switch:
        import eazy.igm
        igm = eazy.igm.Inoue14()
        igm_factor = igm.full_IGM(z, (1+z)*RFull*10**4)



    'DL07'
    DL07=np.zeros((len(sfx),len(dat1[0,:,0]),len(dat1[0,0,:]),len(g)))
    agn_c=np.zeros((len(sfx),len(agn[:,0])))
    GB_T=np.zeros((len(sfx),len(GB[:,0])))

    if extra_bands:
        DL07_ex=np.zeros((len(sfx_extra),len(dat1[0,:,0]),len(dat1[0,0,:]),len(g)))
        agn_c_ex=np.zeros((len(sfx_extra),len(agn[:,0])))
        GB_T_ex=np.zeros((len(sfx_extra),len(GB[:,0])))


    for i,obj in enumerate(dat1[0,:,0]):
        for j,obj in enumerate(dat1[0,0,:]):
            for k,obj in enumerate(g):
                SED=DL07_0[:,i,j,k]*FL*igm_factor
                for p,obj in enumerate(sfx):
                    DL07[p,i,j,k]=convolver(RFullred,SED,FILTERS[p][0],FILTERS[p][1],sfx[p])
                    #DL07[p,i,j,k]=integrate_filter(RFullred,SED,FILTERS[p][0],FILTERS[p][1])




                if extra_bands:
                    for p,obj in enumerate(sfx_extra):
                        DL07_ex[p,i,j,k]=convolver(RFullred,SED,filt_extra[p][0],filt_extra[p][1],sfx_extra[p])
                        #DL07_ex[p,i,j,k]=integrate_filter(RFullred,SED,filt_extra[p][0],filt_extra[p][1],)

    'AGN'

    for i,obj in enumerate(agn[:,0]):
        for p,obj in enumerate(sfx):
            agn_c[p,i]=convolver(RFullred,agn[i,:]*igm_factor,FILTERS[p][0],FILTERS[p][1],sfx[p])
            #agn_c[p,i]=integrate_filter(RFullred,agn[i,:]*igm_factor,FILTERS[p][0],FILTERS[p][1],sfx[p])

        if extra_bands:
            for p,obj in enumerate(sfx_extra):
                agn_c_ex[p,i]=convolver(RFullred,agn[i,:]*igm_factor,filt_extra[p][0],filt_extra[p][1],sfx_extra[p])
    'Optical'

    GB_T=np.zeros((len(sfx),len(GB[:,0])))
    for i,obj in enumerate(GB[:,0]):
        for p,obj in enumerate(sfx):
            #GB_T[p,i]=convolver(RFullred,GB[i,:]*igm_factor,FILTERS[p][0],FILTERS[p][1],sfx[p])
            GB_T[p,i]=integrate_filter(RFullred,GB[i,:]*igm_factor,FILTERS[p][0],FILTERS[p][1])

        if extra_bands:
            for p,obj in enumerate(sfx_extra):
                #GB_T_ex[p,i]=convolver(RFullred,GB[i,:]*igm_factor,filt_extra[p][0],filt_extra[p][1],sfx_extra[p])
                GB_T_ex[p,i]=integrate_filter(RFullred,GB[i,:]*igm_factor,filt_extra[p][0],filt_extra[p][1])



    if extra_bands:
        DL07=np.concatenate([DL07,DL07_ex],axis=0)
        agn_c=np.concatenate([agn_c,agn_c_ex],axis=0)
        GB_T=np.concatenate([GB_T,GB_T_ex],axis=0)

    time_synphot=time.time()-timer_synphot

    if verbose>1:
        print('Synthetic Photometry - ',time_synphot,'s')


    #==========================================================================

    '''
    (Re)Normalising Templates to perform the fit
    '''
    b=10**10*(const.M_sun/const.m_p)*dust_switch
    b1=1.#/np.mean(GB)*stellar_switch
    b2=1./np.mean(agn)*agn_switch


    #-------------------------------------
    #Masking template entries with no photometry
    DL07=np.delete(DL07, mask,axis=0)
    GB_T=np.delete(GB_T,mask,axis=0)
    agn_c=np.delete(agn_c,mask,axis=0)



    sort=np.argsort(wav_ar)
    flux=flux[sort]
    flux_e=flux_e[sort]
    flux_e_orig=flux_e_orig[sort]
    wav_ar=wav_ar[sort]
    DL07=DL07[sort]
    GB_T=GB_T[sort]
    agn_c=agn_c[sort]




    #plt.plot(wav_ar,flux,'ro')
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.show()
    #sys.exit()

    #-------------------------------------
    A=[]
    try:
        CHI=np.zeros_like(DL07[0,:,:,:])
    except:
        return None
    SOL=np.zeros((15,len(DL07[0,:,0,0]),len(DL07[0,0,:,0]),len(DL07[0,0,0,:])))

    for i,obj in enumerate(GB_T[0,:]):
        A.append(b1*GB_T[:,i])


    A.append(b2*agn_c[:,0])
    A.append(b2*agn_c[:,1])
    '''Can insert DL linear combination here'''
    A.append(0*agn_c[:,1]) #Adding a fake row as a placeholder for DL07

    #---------------------------------------------------------------------------
    minsol=np.array([0,0,0])

    def stellar(*params): #Stellar template
        params=np.array(params)
        optt=params*np.array(GB.T)
        optt=np.sum(optt,axis=1)
        return b1*optt


    def agnpl(*params):  #AGN template
        params=np.array(params)
        a=params*b2*agn[1:3,:].T
        a=np.sum(a,axis=1)
        return a


    def full_ir(params,template): #IR template
        params=np.array(params)
        r=[]
        r.append(S(template[0],template[1],template[2]))
        r=np.array(r)
        r=params*r.T
        return b*FL*(np.sum(r,axis=1))


    def full_template(params,template): #All combined
        params=np.array(params)
        s=stellar(params[:steltemp])
        i=full_ir(params[total-irtemp:],template)
        a=agnpl(params[steltemp:steltemp+agntemp])
        return s+i+a

    def full_template_ir(params,template): #All combined
        params=np.array(params)
        s=stellar(params[:steltemp])
        i=full_ir(params[total-irtemp:],template)
        a=agnpl(params[steltemp:steltemp+agntemp])
        return i

    '''
    This computes LIR, given the template function and the solution vector
    '''
    def LumIR(f,param,template=minsol):
        template=np.array(template)
        flux=f(param,template)*10**-26
        luminosity_v=((4*pi*lum_dist**2*flux*(1+z)**-1)*u.erg).si
        wavrange=[find_nearest(RFull,8),find_nearest(RFull,1000)+1]
        luminosity_v=luminosity_v[wavrange[0]:wavrange[1]]
        wavelength=((RFull[wavrange[0]:wavrange[1]]*(1.+0))*u.um).si
        Luminosity_IR=integrate.trapz((luminosity_v*c*wavelength**-2).value,wavelength.value)*(3.839*10**26)**-1
        return Luminosity_IR
    #---------------------------------------------------------------------------
    fit_start=time.time()

    Lir_chi=[]
    Lir_draine_chi=[]
    Mdust_chi=[]
    gamma_chi=[]
    chi2_cov=[]
    fpdr_chi=[]

    ir_cntr=len(wav_ar[wav_ar>3.5])

    chi2_opt_cov=[]

    for i1,obj in enumerate(dat1[0,:,0]):
        for i2,obj in enumerate(dat1[0,0,:]):
            for i3,obj in enumerate(g):
                B=np.copy(A)
                B[-1,:]=b*DL07[:,i1,i2,i3]
                B=np.array(B)
                C=np.copy(B)
                B=B/flux_e
                B=B.T
                try:
                    nnsol, ier = sp.optimize.nnls(B,flux/flux_e)
                    SOL[:,i1,i2,i3]=nnsol
                    CHI[i1,i2,i3]=ier**2
                    if save_covar:
                        Lir_chi.append(LumIR(full_template,nnsol,template=[i1,i2,i3])) #Total LIR
                        Lir_draine_chi.append(LumIR(full_template_ir,nnsol,template=[i1,i2,i3])) #LIR only for the DL07 template
                        #Mdust_chi.append(DG_ratio*nnsol[-1]*b*(const.m_p/const.M_sun)) #Mdust
                        gamma_chi.append(g[i3])
                        fpdr_chi.append(f_pdr(g[i3],Umin[i1]))
                        chi2_cov.append(ier**2)
                        chi2_ir=np.sum(  ((np.sum((nnsol*C.T),axis=1)-flux)**2*(flux_e**2)**-1)[wav_ar>3.5]           )
                        chi2_opt=np.sum(  ((np.sum((nnsol*C.T),axis=1)-flux)**2*(flux_e**2)**-1)[wav_ar<3.5]           )
                        chi2_opt_cov.append(chi2_opt)
                        print('IR',chi2_ir)
                        print('OPT',chi2_opt)
                        print('TOT',ier**2)


                except:
                    nnsol = np.array([0.]*total)

                    ier=-99


    #Find the indices of the best fit solution

    minsol=np.unravel_index(np.argmin(CHI),CHI.shape)
    nnsol=SOL[:,minsol[0],minsol[1],minsol[2]]
    chi2=CHI[minsol[0],minsol[1],minsol[2]]


    #Create matrix with best-fit template
    A[-1]=b*DL07[:,minsol[0],minsol[1],minsol[2]]
    A=np.array(A)
    bestfit=np.copy(A)
    A=A/flux_e
    A=A.T


    SOL_f=SOL.reshape(SOL.shape[0],-1)
    CHI_f=CHI.flatten()


    nfree=1
    #deltaCHI=(CHI_f/(len(flux)-nfree))-(chi2/(len(flux)-nfree))
    deltaCHI=(CHI_f-chi2)
    Delta_chisq_of_interest=11.498779

    CHImask=deltaCHI<Delta_chisq_of_interest

    SOL_13=SOL_f[:,CHImask]

    covmask=np.argwhere(nnsol==0)
    covmask2=nnsol>0
    nnsol_mask=np.delete(nnsol,covmask,None)
    C=np.delete(A,covmask,axis=1)


    Mdust_chi=(DG_ratio*SOL_f[-1,:]*b*(const.m_p/const.M_sun)) #Mdust

    #------------------------------------------------------------------
    #Functions
    '''
    def stellar(*params): #Stellar template
        params=np.array(params)
        optt=params*np.array(GB.T)
        optt=np.sum(optt,axis=1)
        return b1*optt


    def agnpl(*params):  #AGN template
        params=np.array(params)
        a=params*b2*agn.T
        a=np.sum(a,axis=1)
        return a


    def full_ir(*params): #IR template
        params=np.array(params)
        r=[]
        r.append(S(minsol[0],minsol[1],minsol[2]))
        r=np.array(r)
        r=params*r.T
        return b*FL*(np.sum(r,axis=1))


    def full_template(*params): #All combined
        params=np.array(params)
        s=stellar(params[:steltemp])
        i=full_ir(params[total-irtemp:])
        a=agnpl(params[steltemp:steltemp+agntemp])
        return s+a+i


    def LumIR(f,param):
        flux=f(*param)*10**-26
        luminosity_v=((4*pi*lum_dist**2*flux*(1.+z)**-1)*u.erg).si
        wavrange=[find_nearest(RFull,8),find_nearest(RFull,1000)+1]
        luminosity_v=luminosity_v[wavrange[0]:wavrange[1]]
        wavelength=((RFull[wavrange[0]:wavrange[1]]*(1.+0))*u.um).si
        Luminosity_IR=integrate.trapz((luminosity_v*c*wavelength**-2).value,wavelength.value)*(3.839*10**26)**-1

        return Luminosity_IR
    '''
    #---------------------------------------------------
    chi2_reduced=chi2/(len(flux)-nfree)
    chi2_manual=np.sum((np.sum((nnsol*bestfit.T),axis=1)-flux)**2*(flux_e**2)**-1)
    #---------------------------------------------------
    #Finding initial values for quantities of interest
    #Lir_draine=LumIR(full_ir,nnsol[(total-irtemp):])
    #Lir_total=LumIR(full_template,nnsol)
    #Lir_stellar=LumIR(stellar,nnsol[:steltemp])

    Lir_draine=LumIR(full_ir,nnsol[(total-irtemp):],template=minsol)
    Lir_total=LumIR(full_template,nnsol,template=minsol)
    #Lir_stellar=LumIR(stellar,nnsol[:steltemp])



    Lagn=(Lir_total-Lir_draine)
    fagn=(Lir_total-Lir_draine)/Lir_total

    if np.isfinite(fagn):
        fagn=fagn
    else:
        fagn=0.
    SFR=Lir_draine*10**-10


    Mgas_nnls=nnsol[-1]*b*(const.m_p/const.M_sun)

    Mdust=DG_ratio*Mgas_nnls




    #---------------------------------------------------------------------------

    dl=cosmo.luminosity_distance(z).to(u.cm)
    conv=((1*u.mJy).to(u.erg/u.s/u.cm**2/u.Hz))*4*np.pi*dl**2/(1+z)

    fnu_units = u.erg/u.s/u.cm**2/u.Hz
    mJy_to_cgs = u.mJy.to(u.erg/u.s/u.cm**2/u.Hz)
    fnu_scl = 10**(-0.4*(ABZP-23.9))*1*mJy_to_cgs

    template_fnu_units=(1*u.solLum / u.Hz)
    to_physical = fnu_scl*fnu_units*4*np.pi*dl**2/(1+z)
    to_physical /= (1*template_fnu_units).to(u.erg/u.second/u.Hz)


    coeffs_rest = (b1**-1*nnsol[:12].T*to_physical).T

    # Remove unit (which should be null)
    coeffs_rest = np.array(coeffs_rest)

    mass = coeffs_rest.dot(t_param['mass'])
    sfr = coeffs_rest.dot(t_param['sfr'])
    #print(np.log10(mass))
    print('cat_mass',np.log10(Mstar))


    #Compute stellar mass from the K-band, using L/M=0.6

    idx_k=find_nearest(RFull,2.2)
    L_K=(stellar(*nnsol[:(steltemp)])*igm_factor*conv)[idx_k]*(const.c/(2.2*u.um))
    L_K=((L_K).to(u.erg/u.s)).to(u.solLum).value
    mass_K=L_K*0.6
    #print('mass K',np.log10(mass_K))

    #---------------------------------------------------------------------------

    U=Lir_draine/(125*Mdust)

    reddest_flux=DL07[-1,:,:,:].flatten()
    reddest_flux/=np.median(reddest_flux)
    reddest_flux=reddest_flux[deltaCHI<Delta_chisq_of_interest]

    red_sigma=Mdust*((10**np.std(np.log10(reddest_flux)))-1)





    #Metallicity from Manucci+10 eq5
    #------------------------------------------------------
    if Mstar==-99:
        Mgas=deltaGDR=-99
    else:
        if use_own_stellar_mass:
            Mgas,deltaGDR=mgas_metallicity_manucci(mass,Mdust,SFR)
        else:
            Mgas,deltaGDR=mgas_metallicity_manucci(Mstar,Mdust,SFR)
    #------------------------------------------------------

    cov_start=time.time()

    Mdust_cov=[]
    Lir_total_cov=[]
    Lir_draine_cov=[]
    eLir_total=eLir_draine=eLagn=eMD=eMG=efagn=-99
    max_allowed = 5
    attempt = 0
    timer=1
    loop_time_start=time.time()
    while (eLir_total<1.) or np.isnan(eLir_total) or (eMD<1.) or np.isnan(eMD):
        attempt+=1
        if attempt>max_allowed or time.time()-loop_time_start>timer:
            print('multivar_gauss: max attempts exceeded: TIMEOUT ERROR')
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
            break

        try:
            nnsol_cov_sum=np.sum(nnsol_cov[:,total-irtemp:],axis=1)
            Mdust_cov=DG_ratio*nnsol_cov_sum*b*(const.m_p/const.M_sun)
            coeffs_rest = (b1**-1*nnsol_cov[:,:12].T*to_physical).T
            coeffs_rest = np.array(coeffs_rest)
            #mass_cov = coeffs_rest.dot(t_param['mass'])
            mass_cov=np.tensordot(coeffs_rest,t_param['mass'],axes=1)

            for x,obj in enumerate(nnsol_cov[:,0]):
                Lir_total_cov.append(LumIR(full_template,nnsol_cov[x,:],template=minsol))
                Lir_draine_cov.append(LumIR(full_ir,nnsol_cov[x,(total-irtemp):],template=minsol))

            Lagn_cov=np.array(Lir_total_cov)-np.array(Lir_draine_cov)
            fagn_cov=np.array(Lagn_cov)/np.array(Lir_total_cov)

            eLagn=mean_confidence_interval_fast_error(Lagn_cov,Lagn)
            eLir_draine=mean_confidence_interval_fast_error(Lir_draine_cov,Lir_draine)
            eLir_total=mean_confidence_interval_fast_error(Lir_total_cov,Lir_total)
            eMD=mean_confidence_interval_fast_error(Mdust_cov,Mdust)
            efagn=mean_confidence_interval_fast_error(fagn_cov,fagn)
            eMG=eMD*deltaGDR
            eMstar=mean_confidence_interval_fast_error(mass_cov,mass)

        except Exception as e:
            print(e)
            print('multivar_gauss: fail')
            break

    #----------------------------------------------------------------------------

    try:
        print(doesntexist)
        Lir_chi=[]
        Lir_draine_chi=[]
        Mdust_chi=[]
        gamma_chi=[]
        for x,obj in enumerate(SOL_13[0,:]):
            Lir_chi.append(LumIR(full_template,SOL_13[:,x]))
            Lir_draine_chi.append(LumIR(full_ir,SOL_13[(total-irtemp):,x]))
            Mdust_chi.append(DG_ratio*SOL_13[-1,x]*b*(const.m_p/const.M_sun))

        Lir_chi=np.array(Lir_chi)
        Mdust_chi=np.array(Mdust_chi)
        Lir_total_cov=np.array(Lir_total_cov)
        Mdust_cov=np.array(Mdust_cov)
        CIL=mean_confidence_interval_fast(np.log10(Lir_chi))
        CIM=mean_confidence_interval_fast(np.log10(Mdust_chi))


        Lir_med=10**CIL[0]
        Mdust_med=10**CIM[0]

        Lir_med_err=np.mean([10**CIL[0]-10**CIL[1],10**CIL[2]-10**CIL[0]])
        Mdust_med_err=np.mean([10**CIM[0]-10**CIM[1],10**CIM[2]-10**CIM[0]])
    except:
        Lir_med=-99
        Mdust_med=-99
        Lir_med_err=-99
        Mdust_med_err=-99


    sU=(1/125)*np.sqrt((eLir_draine/Mdust)**2+(Lir_draine*eMD/Mdust**2)**2)


    CHI_13=CHI_f[CHImask]

    if eMD>0:
        eMD=(eMD**2+red_sigma**2)**0.5
    else:
        eMD=red_sigma

    emd2=np.std(np.array(Mdust_chi)[deltaCHI<Delta_chisq_of_interest])


    if save_covar:

        cov_data=np.array([chi2_cov,Lir_chi,Lir_draine_chi,Mdust_chi,gamma_chi,fpdr_chi])
        cov_names=['chi2','Lir_total','Lir_draine','MD','gamma','fpdr']
        covar=Table(cov_data.T,names=cov_names,dtype=len(cov_names)*[float])
        covar.write(covarloc+str(int(P[galaxy_index,0])) + '_covar_all_chi.fits',overwrite=True)

    if verbose:
        print('------------------------------------------')
        print ('ID:',int(P[galaxy_index,0]))
        print ('LIR:',"%.4g" % (Lir_total),'+-',"%.4g" % (eLir_total),'Lsol')
        print ('Mdust:',"%.4g" % Mdust ,'+-',"%.4g" % eMD,'Msol' )
        print ('Mstar:',"%.4g" % mass ,'+-',"%.4g" % eMstar,'Msol' )
        print ('chi2_nu:',"%.4g" % chi2_reduced)
        print('------------------------------------------')

    #--------------------------------------------------
    def rad_flux(lam,method='delv20',alpha=-0.75,use_K=False):
        '''
        Calculate the radio contribution:
        Two methods available:
        Delvecchio+20, with Lir and Mstar dependence - delv20
        Delhaize+17, only with Lir dependence - delhz17
        If no stellar mass available, can use the K-band prediction instead - use_K
        '''
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
    #--------------------------------------------------

    if save_sed:

        if radio:
            RADIO=rad_flux(RFull*(1.+z),method=radio_method,alpha=-0.75,use_K=False)
        else:
            RADIO=np.zeros_like(RFull)

        toSED=np.array([RFull*(1.+z),stellar(*nnsol[:(steltemp)])*igm_factor,agnpl(nnsol[steltemp:steltemp+agntemp])*igm_factor,full_ir(nnsol[(total-irtemp):],minsol)*igm_factor,full_template(nnsol,minsol)*igm_factor,RADIO])
        tableSED=Table(toSED.T,names=['lambda','stellar','AGN','IR','Total','Radio'],dtype=[float,float,float,float,float,float])
        tableSED.write(sedloc+str(int(P[galaxy_index,0])) + ".fits",overwrite=True)

    #--------------------------------------------------

    if not multithread:
        fig= plt.figure(figsize=(10,5))
        ax = fig.add_subplot(1,1,1)
        textypos=0.9
        textxpos=0.6
        textsep=0.08



        #for i,_ in enumerate(nnsol_cov[:,0]):
        #    plt.plot(RFull*(1.+z),full_template(*nnsol_cov[i,:]),'gray',lw=3,alpha=0.1,zorder=-1)
        #plt.plot(RFull*(1.+z),full_template(*nnsol),'k',label='Total',lw=3,alpha=0.8)

        SF=stellar(*nnsol[:(steltemp)])*igm_factor
        AGN=agnpl(nnsol[steltemp:steltemp+agntemp])*igm_factor
        IR=full_ir(nnsol[(total-irtemp):],minsol)*igm_factor

        if radio:
            #RADIO=rad_flux(RFull*(1.+z))
            RADIO=rad_flux(RFull*(1.+z),method=radio_method,alpha=-0.75,use_K=False)
            TOTAL=SF+AGN+IR+RADIO
        else:
            TOTAL=SF+AGN+IR

        #radio_bands=np.array([30000,200000])
        #radio_points=np.array([1.1, 22.4])*10**-3
        #e_radio_points=np.array([1.1, 6.4])*10**-3

        plt.fill_between(RFull*(1.+z),0,SF,color='royalblue',alpha=0.2,label='Stellar')
        plt.fill_between(RFull*(1.+z),0,AGN,color='g',alpha=0.2,label='AGN')
        plt.fill_between(RFull*(1.+z),0,IR,color='maroon',alpha=0.2,label='Dust')


        plt.plot(RFull*(1.+z),TOTAL,'k',label='Total',lw=3,alpha=0.6,zorder=10)


        #plt.plot(wav_ar,np.sum(nnsol*bestfit.T,axis=1),'bo')



        points=((flux/flux_e_orig)>=3)


        plt.errorbar(wav_ar[points],flux[points],yerr=flux_e_orig[points],color='red',fmt='s',capsize=5,capthick=1,ms=12,markerfacecolor='white',mew=2,barsabove=True)
        plt.scatter(wav_ar[~points],(flux+3*flux_e_orig)[~points],marker=r'$\downarrow$',s=300,color='red',zorder=11)

        #if radio:
        #   points_radio=((radio_points/e_radio_points)>=3)
        #   plt.errorbar(radio_bands[points_radio],radio_points[points_radio],yerr=e_radio_points[points_radio],color='blue',fmt='o',capsize=5,capthick=1,ms=12,markerfacecolor='white',mew=2,barsabove=True)
        #   plt.scatter(radio_bands[~points_radio],(radio_points+3*e_radio_points)[~points_radio],marker=r'$\downarrow$',s=300,color='b',zorder=11)


        plt.text(0.02, 0.85,'ID '+str(int(P[galaxy_index,0])), color='k',fontsize=20,transform=ax.transAxes)
        plt.text(0.02, 0.75,r'z={:.2f}'.format(z), color='k',fontsize=20,transform=ax.transAxes)
        plt.text(0.02, 0.65,r'logf={:.2f}'.format(np.log10(Mgas/Mstar)), color='k',fontsize=20,transform=ax.transAxes)

        ylabel(r'$f_{\nu}$ [mJy]',fontsize=25)
        xlabel(r'$\lambda_{obs}$ $[\mu m]$',fontsize=25)
        legend(fontsize=12)
        ylim(10**-4,10**3)
        xlim(.5,10**5.7)
        grid(alpha=0.4)
        yscale('log')
        xscale('log')
        plt.tight_layout()

        if diagplot:
            show()

    #----------------------------------------------------------------------------------------

    if verbose:
        print('------------------------------------------')
        print ('Finished fitting ID:',P[galaxy_index,0])
        print('------------------------------------------')

    R=np.array([int(P[galaxy_index,0]),Lir_total,eLir_total,Mdust,eMD,z,chi2_reduced,fagn,efagn,
    lastdet,Mgas,eMG,deltaGDR,attempt,Mstar,100*Mdust/Mstar,
    Mgas/Mstar,Lir_med,Lir_med_err,Mdust_med,Mdust_med_err,Umin[minsol[0]],minsol[1],
    g[minsol[2]],U,sU,Lagn,eLagn,Lir_draine,eLir_draine,mass,eMstar,mass_K,chi2,nfilt])

    nnsol=np.array(nnsol)



    if save_table:
        att=0
        if multithread:
            with lock:
                t0=Table.read(table_out)
                t0.add_row(R)
                t0.write(table_out,overwrite=True)
        else:
            t0=Table.read(table_out)
            t0.add_row(R)
            t0.write(table_out,overwrite=True)
    return None


objects=range(len(G[:,0]))

if not multithread:
    otp=[]
    for i,obj in enumerate(objects):
        otp.append(galproc(i))


if multithread:
    from multiprocessing import Pool
    from multiprocessing import Process, Lock
    import multiprocessing
    print('Begin multithreading')
    lock = Lock()
    pool = Pool(multiprocessing.cpu_count())
    print(multiprocessing.cpu_count(),'threads utilised')
    mp_out= np.array(pool.map(galproc,objects))
    mp_out=np.array(mp_out)
    print ('Code took',time.time()-timer_start, 's to run')


if save_table:
    table_final=Table.read(table_out)
    table_final=table_final[1:] #Remove placeholder row
    table_final.write(table_out,overwrite=True)

if save_fig:
    from plotter import *

sys.exit()
