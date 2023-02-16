import os
import sys
import time
from os.path import exists

import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM

from . import utils
from stardust.filter_importer import *

loc = os.path.dirname(os.path.realpath(__file__))
templ_loc = os.path.join(loc,'templates')

try:
    from tqdm import tqdm
    HAS_TQDM = True
except:
    HAS_TQDM = False

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

c=const.c
Lsol=(const.L_sun.cgs).value
Msol=(const.M_sun).value
pi=np.pi


steltemp=12 #Default is 12
agntemp=2 #Default is 2
irtemp=1 #Default is 1

total=steltemp+agntemp+irtemp

#Some legacy params
detection_threshold=3 #In sigma
ir_detections=3 
ir_cutoff=20

__version__ = "1.0.2"

if (sys.version_info <= (3,0)):
    print(f'Python Version {sys.version[:6]}, is unsupported. Please use Python 3.9+')
    sys.exit()
    
if (sys.version_info >= (3, 5)) and (sys.version_info <= (3, 7)):
    old_ver = True
    print(f'Detected Python {sys.version[:6]}, switching multiprocessing packages')
    
if (sys.version_info > (3, 7)):
    old_ver = False


#===============================================================================

class ctf(object):

    def __init__(self, idx = None,config_file=None,zeropoint_file=None,**kwargs):
       
        self.kwargs = kwargs

        try:
            if kwargs['custom_umin_indices'] is None:
                pass
        except:
            kwargs = {'custom_gamma': None,
                        'custom_qpah_indices':None,
                            'custom_umin_indices': None}


        self.config_file = config_file
        self.zp_file = zeropoint_file

        self.filters_all = filters

        self.read_config()
        self.read_catalogue(idx=idx)
        self.make_template_grid(**kwargs)
        self.prepare_products()

        return None

    def read_config(self):
        """
        READ THE CONFIG FILE
        """
        self.config = {}

        f = open(self.config_file,'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            if not line.strip().startswith('#'):
                lsplit = line.split()
                if len(lsplit)>2:
                    sys.exit('Incorrect formatting of the config file, check example')
                #if lsplit[1].isnumeric():
                if lsplit[1].replace('.','',1).isnumeric():
                    lsplit[1]=float(lsplit[1])
                self.config[lsplit[0]]=lsplit[1]
        print(f'READ CONFIG FILE {self.config_file}')

        return None

    def read_catalogue(self, idx = None, cleanup=True):
        """
        READ THE INPUT CATALOGUE
        ADD THE NECESSARY FILTERS TO FILTER LIST
        """

        #Function too long, maybe break up into a few

        try:
            self.data = Table.read(self.config['CATALOGUE'])
        except:
            self.data = Table.read(self.config['CATALOGUE'], format='ascii')
        
        if idx is not None:
            self.data = self.data[idx[0]:idx[1]]
            
        print(f'Read a catalogue with {len(self.data)} objects')

        bnds=np.loadtxt(self.config['BANDS_FILE'],dtype='str')
        prms=np.loadtxt(self.config['PARAM_FILE'],dtype='str')

        if self.config['VERBOSE']>0:
            print('Read band file: ' + self.config['BANDS_FILE'])
            print('Read param file: ' + self.config['PARAM_FILE'])

        if len(bnds[0,:])>2:
            print('Using Legacy Format')
            self.f_numbers=bnds[:,0].tolist() #Band Filter ID
            band_names=bnds[:,1].tolist() #Band Name
            err_band_names=bnds[:,2].tolist() #Error
        
        else:
            print('Using eazy-py Format')
            '''Load in the file as dictionary to remove potential duplicates'''
            filter_dict = {}
            lines = open(self.config['BANDS_FILE']).readlines()

            for line in lines:
                if line.startswith('#'):
                    continue
                spl = line.split()
                key = spl[0]
                filter_dict[key] = spl[1]

            self.f_numbers = []
            band_names = []
            err_band_names = []
            for i,val in enumerate(filter_dict.values()):

                fname = val[1:]
                name = list(filter_dict.keys())[i]
                if fname not in self.f_numbers:
                    self.f_numbers.append(fname)

                if val.startswith('F'):
                        band_names.append(name)

                if val.startswith('E'):
                    err_band_names.append(name)
                  
        param_names=prms[:,1].tolist() #Get list of parameters (z,Mstar etc.)

        if param_names[2]=='None':
            self.data['Mstar']=-99 #If catalogue does not have its own input Mstar 
            param_names[2]='Mstar'


        filt_id=np.int_(self.f_numbers).tolist()


        #Append the filters to be used later on

        self.filters=[]
        self.sfx=[]
        self.fnames=[]

        for id in filt_id:
            fltr=np.array([self.filters_all[id-1].wave*10**-4,self.filters_all[id-1].throughput])
            self.filters.append(fltr)
            lspl=filters[id-1].name.split()
            self.fnames.append(lspl[0])
            for index,s in enumerate(lspl):
                if 'lambda' in s:
                    wav_str = lspl[index+1]
                    #if 'um' in wav_str:
                    #    wav_str = wav_str[:-2]
                    wav_str = wav_str.strip('um')
                    #wav_num =  ''.join(filter(str.isdigit, wav_str) ) #Make sure that only numbers are used
                    self.sfx.append(wav_str)
                    #self.sfx.append(lspl[index+1])
        self.sfx=1e-4*np.float_(self.sfx) #Convert central wavelength to microns


        assert len(self.fnames)==len(band_names)==len(err_band_names),'ERROR: Filter range does not match available photometry, check band names'

        if self.config['VERBOSE']>0:
            if self.config['VERBOSE']>1:
                print(f'Detected the following filters...')
                print('--------------------------------')
                for name in self.fnames:
                    print(name)
                    print('--------------------------------')

        print(f'Detected the following parameters {param_names}')
        print('Input Bands: ok')
        print('Input Filters: ok')


        self.fnu = np.array((np.array(self.data[band_names])).tolist())
        self.efnu = np.array((np.array(self.data[err_band_names])).tolist())
        self.wav = np.tile(self.sfx,(len(self.fnu[:,0]),1))

        
        if self.zp_file is not None:
            print(self.zp_file)
            self.zp = self.apply_zeropoint()
            if len(self.zp)==len(self.fnu[0,:]):
                print('Applying Zeropoint corrections')
                self.fnu*=self.zp
                self.efnu*=self.zp
            else:
                print('Length mismatch in the zeropoints and bands in the catalogue, check the zeropoint file')


        #self.param = np.array((np.array(self.data[param_names])).tolist())

        self.id = self.data[param_names[0]]
        self.zcat = self.data[param_names[1]]
        self.mstar_cat = self.data[param_names[2]]

        self.lnu_to_fnu = (1.+self.zcat)*(4.*pi*(cosmo.luminosity_distance(self.zcat)).to(u.cm).value**2)**-1

        if self.config['EXTRA_BANDS']:
            self.add_extra_bands()
            self.fnu=np.concatenate((self.fnu,self.fnu_extra),axis=1)
            self.efnu=np.concatenate((self.efnu,self.efnu_extra),axis=1)
            self.wav=np.concatenate((self.wav,self.wav_extra),axis=1)


        self.fnu*=self.convert_to_mjy(self.config['FLUX_UNIT'])
        self.efnu*=self.convert_to_mjy(self.config['FLUX_UNIT'])

        self.NOBJ = len(self.fnu)

        coord_names = np.array([['RA', 'ra', 'ALPHA_J2000', 'alpha'] ,
                                [ 'DEC', 'dec', 'DELTA_J2000', 'delta']])

        for col in self.data.colnames:
            if col in coord_names[0]:
                ra_col = col
            if col in coord_names[1]:
                dec_col = col
        
        try:
            self.ra = self.data[ra_col]
            self.dec = self.data[dec_col]
        except:
            print('Could not automatically assign the ra,dec columns')
            self.ra = np.zeros_like(len(self.data))
            self.dec = np.zeros_like(len(self.data))

        print('Data Loaded Successfully')
        print('Objects to Fit:',len(self.fnu[:,0]))

        return None


    def convert_to_mjy(self,unit):
        """
        Convert Input Units to mJy
        """

        units={'Jy':(u.Jy).to(u.mJy),
                'mJy':1.,
                'uJy':(u.uJy).to(u.mJy),
                'nJy':(u.nJy).to(u.mJy)}
        try:
            conv=units[unit]
        except:
            sys.exit(f'Unrecognised unit - {unit}, accepted formats: Jy, uJy, mJy, nJy')
        return conv

    def add_extra_bands(self):
        """
        ADD THE CUSTOM FILTERS AND BANDS TO CATALOGUE
        """

        if self.config['EXTRA_BANDS']: #Redundant but just to be safe
            xtra_bnds=np.loadtxt(self.config['EXTRA_BANDS_FILE'],dtype='str')
        try:
            wavelength_extra=xtra_bnds[:,0].tolist()
            bands_extra=xtra_bnds[:,1].tolist()
            err_bands_extra=xtra_bnds[:,2].tolist()
        except:
            wavelength_extra=[np.array(xtra_bnds[0]).tolist()]
            bands_extra=[np.array(xtra_bnds[1]).tolist()]
            err_bands_extra=[np.array(xtra_bnds[2]).tolist()]     
        if self.config['VERBOSE']==1:
            print(f'Added extra bands from ' + self.config['EXTRA_BANDS_FILE'])

        self.fnu_extra=np.array((np.array(self.data[bands_extra])).tolist())
        self.efnu_extra=np.array((np.array(self.data[err_bands_extra])).tolist())
        self.wav_extra=np.array((np.array(self.data[wavelength_extra])).tolist())
    
        return None
        

    def filter_diagnostic(self):
        """
        Plots all of the filter curves that will be used in the fitting
        """

        import math
        xlen=6
        ylen=math.ceil(len(self.filters)/6)
        fig=plt.figure(figsize=(3*xlen,3*ylen))

        for f,_ in enumerate(self.filters[:]):
            ax = fig.add_subplot(ylen,xlen,f+1)
            ax.plot(self.filters[f][0],self.filters[f][1])
            ax.axvline(self.sfx[f],color='r')
            ax.text(0.2,0.2,np.round(self.sfx[f],2),transform=ax.transAxes)

        plt.show()

        return None

    def prepare_products(self):
        '''
        Create placeholder output tables and folders
        '''

        self.usephot = np.zeros((self.NOBJ,self.fnu.shape[1]))
        self.best_coeffs = np.zeros((self.NOBJ,total))
        self.best_ir_idx = np.zeros((self.NOBJ,3))

        self.lir = np.zeros(self.NOBJ)
        self.lir_agn = np.zeros(self.NOBJ)
        self.lir_tot = np.zeros(self.NOBJ)
        self.sfr = np.zeros(self.NOBJ)

        self.fagn = np.zeros(self.NOBJ)
        self.e_fagn = np.zeros(self.NOBJ)
        self.mdust = np.zeros(self.NOBJ)
        self.deltaGDR = np.zeros(self.NOBJ)
        self.mgas = np.zeros(self.NOBJ)
        self.mstar = np.zeros(self.NOBJ)
        self.sfr_optical = np.zeros(self.NOBJ)
        self.av = np.zeros(self.NOBJ)

        self.e_lir_tot = np.zeros(self.NOBJ)
        self.e_lir_agn = np.zeros(self.NOBJ)
        self.e_lir = np.zeros(self.NOBJ)
        self.e_sfr = np.zeros(self.NOBJ)
        self.e_mdust = np.zeros(self.NOBJ)
        self.e_mgas = np.zeros(self.NOBJ)
        self.e_mstar = np.zeros(self.NOBJ)
        self.e_sfr_optical = np.zeros(self.NOBJ)
        self.e_av = np.zeros(self.NOBJ)

        self.chi2 = np.zeros(self.NOBJ)
        self.lastdet = np.zeros(self.NOBJ)
        self.nfilt = np.zeros(self.NOBJ)
        
        #DL07 params
        self.avu = np.zeros(self.NOBJ)
        self.e_avu = np.zeros(self.NOBJ)
        self.umin = np.zeros(self.NOBJ)
        self.qpah = np.zeros(self.NOBJ)
        self.gamma = np.zeros(self.NOBJ)

        for i in range(self.NOBJ):
            sn = self.fnu[i]/self.efnu[i]
            rf = (1+self.zcat[i])*self.sfx
            try: 
                self.lastdet[i] = (rf[sn>=3][-1])
            except:
                self.lastdet[i] = - 99.
            
        
        
        if self.config['SAVE_TABLE']:
            try:
                os.makedirs(self.config['PATH'])
            except Exception as e:
                print(e)
                pass


            #table_out=self.config['PATH']+self.config['OUTPUT_NAME']+'.fits'
            #self.make_output_table(table_out)

       
        if self.config['SAVE_FIGURE']:
            figloc=self.config['PATH']+'figures/'
            try:
                os.makedirs(figloc)
            except Exception as e:
                print(e)
                pass

        if self.config['SAVE_SED']:
            sedloc=figloc=self.config['PATH']+'seds/'
            try:
                os.makedirs(sedloc)
            except Exception as e:
                print(e)
                pass

        if self.config['SAVE_COVAR']:
            covarloc=figloc=self.config['PATH']+'covar/'
            try:
                os.makedirs(covarloc)
            except Exception as e:
                print(e)
                pass

        return None


    def make_template_grid(self,diagnostic_optical=False,**kwargs):
        '''
         Load optical templates, define common wavelength range
         Load AGN templates
         Load IR templates
        '''

        OT=Table.read(f'{templ_loc}/narrow_no_lines_ext_fix.fits')
        OT_arr = (np.lib.recfunctions.structured_to_unstructured(OT.as_array())).T

        self.templ = np.array(OT['lambda']) #Template wavelength array in microns
        self.tempnu = c.value*(self.templ*1e-6)**-1 #In Hz
        self.flamfnu = (self.templ**2)/c.value #Flambda to Fnu conversion shaped like template wavelength array


        self.optical_param=Table.read(f'{templ_loc}/fsps_QSF_12_v3_narrow.param.fits')


        #ir_range=np.array([find_nearest(self.templ,8),find_nearest(self.templ,1000)+1])

        self.optical_grid = OT_arr[1:,:]

        self.optical_grid = self.optical_grid.T

        if diagnostic_optical:
            fig = plt.figure(figsize=(10,5))
            ax = fig.add_subplot(1,1,1)
            ax.loglog()
            ax.set_ylim(1e-5*np.max(self.optical_grid),10*np.max(self.optical_grid))
            ax.set_xlim(5e-2,1e3)
            ax.set_xlabel(r'$\lambda$ [$\mu$m]')
            ax.set_ylabel(r'$f_\nu$')
            for i,_ in enumerate(self.optical_grid[0,:]):
                ax.plot(self.templ,self.optical_grid[:,i])
            plt.show()


       
        if self.config['QSO']:
            self.optical_grid*=0
            QSO,qso_wave=self.add_qso_templates()

            for i,_ in enumerate(QSO[:,0]):
                qso_func=sp.interpolate.interp1d(qso_wave*10**-4,QSO[i,:],bounds_error=False,fill_value=0,kind='cubic')
                qso=qso_func(self.templ)*self.flamfnu
                self.optical_grid[:,i]=qso

        self.add_agn_templates()

        self.add_ir_templates(**kwargs)

        self.dustnorm = (10**10*(const.M_sun/const.m_p)).value



        return None


    def add_agn_templates(self):

        try:
            self.templ
        except:
            print('Common wavelength grid is not defined, run make_template_grid first')
            return None

        AGN_l=Table.read(f'{templ_loc}/AGN/LoLum_AGN_GM.fits')
        AGN_h=Table.read(f'{templ_loc}/AGN/HiLum_AGN_GM.fits')

        agn_h=sp.interpolate.interp1d(AGN_h['wave']*10**-4,AGN_h['flux'],bounds_error=False,fill_value=0,kind='cubic')
        agn_l=sp.interpolate.interp1d(AGN_l['wave']*10**-4,AGN_l['flux'],bounds_error=False,fill_value=0,kind='cubic')

        agn_high=agn_h(self.templ)*self.flamfnu
        agn_low=agn_l(self.templ)*self.flamfnu

        self.agn_grid=np.zeros((len(agn_low),2))

        self.agn_grid[:,0]=agn_high
        self.agn_grid[:,1]=agn_low


        return None


    def add_ir_templates(self,diagnostic_plot=False,**kwargs):
        """
        CONFIGURE AND LOAD Draine & Li 2014 templates
        """

        try:
            self.templ
        except:
            print('Common wavelength grid is not defined, run make_template_grid first')
            return None

        model1_full = Table.read(f'{templ_loc}/DL07/model1_full.fits')
        model2_full = Table.read(f'{templ_loc}/DL07/model2_full.fits')
        lam = fits.open(f'{templ_loc}/DL07/lambda.fits')

        diffuse=np.array(model1_full['MODEL1'][0]) #DL diffuse part as 3d array
        pdr=np.array(model2_full['MODEL1'][0])  #DL pdr part as 3d array
        dl_lambda=lam[0].data

        self.dg_ratio = 0.01 #Dust-to-gas ratio for the templates, see DL07 and DL14 for more details

        self.templ_umin = np.array([0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,1.000,1.200,
                                    1.700,2.000,3.000,4.000,5.000,7.000,8.000,10.00,12.00,15.00,
                                    20.00,25.00,30.00,35.00,40.00,50.00])


        if not self.config['USE_COLD_DL']:
            custom_umin_indices = range(len(self.templ_umin[self.templ_umin>=0.800]))

        if kwargs['custom_umin_indices'] is None:
            self.umin_indices = range(len(self.templ_umin))
        else:
            self.umin_indices = kwargs['custom_umin_indices']
            self.templ_umin = self.templ_umin[self.umin_indices]

        if kwargs['custom_qpah_indices'] is None:
            self.q_indices = range(len(diffuse[0,0,:]))
        else:
            self.q_indices = kwargs['custom_qpah_indices']

        if kwargs['custom_gamma'] is None:
            g1=np.array([0,0.001,0.0025,0.005,0.0075])
            #g2=np.arange(0.01,0.1,0.01) #Too many gamma values, reduce
            g2=np.arange(0.01,0.1,0.05)
            g3=np.array([0.2,0.35,0.5])
            self.templ_gamma=np.concatenate((g1,g2,g3))

        else:
            self.templ_gamma = kwargs['custom_gamma']

        mesh = np.ix_(range(diffuse.shape[0]),self.umin_indices,self.q_indices)
        diffuse = diffuse[mesh]
        pdr = pdr[mesh]

    
    
        def S(i,j,k):
            """
            Combine diffuse and pdr for parameters of interest
            Interpolate onto the common wavelength grid in log space
            """
            sed=((1.-self.templ_gamma[k])*diffuse[:,i,j]+self.templ_gamma[k]*pdr[:,i,j])
            ss=sp.interpolate.interp1d(np.log10(dl_lambda), np.log10(sed), bounds_error=False,fill_value='extrapolate', kind='linear')
            sed=10**ss(np.log10(self.templ))
            sed=np.array(sed)
            return sed

        self.ir_grid=np.zeros((len(self.templ),len(self.templ_umin),len(self.q_indices),len(self.templ_gamma)))


        for i,_ in enumerate(self.templ_umin):
            for j,_ in enumerate(self.q_indices):
                for k,_ in enumerate(self.templ_gamma):
                    self.ir_grid[:,i,j,k]=S(i,j,k)

        if diagnostic_plot:
            flat_grid = self.ir_grid.reshape(self.ir_grid.shape[0], np.product(self.ir_grid.shape[1:]))
            random_draws = np.random.randint(0,flat_grid.shape[1],100) #Draw 100 random templates
            flat_grid = flat_grid[:,random_draws]
            fig = plt.figure(figsize=(10,5))
            ax = fig.add_subplot(1,1,1)
            ax.loglog()
            ax.set_xlim(1e-1,1e5)
            ax.set_xlabel(r'$\lambda$ [$\mu$m]')
            ax.set_ylabel(r'$f_\nu$')

            for i,_ in enumerate(flat_grid[0,:]):
                ax.plot(self.templ,flat_grid[:,i])

            plt.show()
       
        return None


    def add_qso_templates(self):
        """
        GET SHEN+16 OPTICAL QUASAR TEMPLATES INSTEAD OF BRAMMER+08
        """
        qso_templ1=Table.read(f'{templ_loc}/QSO/shen2016_ext_Av0.0.fits')
        qso_templ2=Table.read(f'{templ_loc}/QSO/shen2016_ext_Av0.1.fits')
        qso_templ3=Table.read(f'{templ_loc}/QSO/shen2016_ext_Av0.2.fits')
        qso_templ4=Table.read(f'{templ_loc}/QSO/shen2016_ext_Av0.3.fits')
        qso_templ5=Table.read(f'{templ_loc}/QSO/shen2016_ext_Av0.4.fits')

        qso_wave=qso_templ1['wave']

        QSO=np.zeros((5,len(qso_wave)))

        QSO[0,:]=qso_templ1['flux']
        QSO[1,:]=qso_templ2['flux']
        QSO[2,:]=qso_templ3['flux']
        QSO[3,:]=qso_templ4['flux']
        QSO[0,:]=qso_templ5['flux']

        return QSO,qso_wave



    def make_extra_filters(self,wavelength,w,dx):
        '''
         Create custom filter arrays
        '''

        lmbd=np.array(wavelength)

        fltrx=[]
        fltry=[]
        for i,_ in enumerate(lmbd):
            c=lmbd[i]
            fltr = np.arange(c-w, c+w, dx)
            fltrtop = np.ones_like(fltr)
            fltrtop[(fltr<c-w) | (fltr>c+w)]=0
            fltrx.append(fltr)
            fltry.append(fltrtop)

        fltrx = np.array(fltrx)
        fltry = np.array(fltry)
        of = np.array([fltrx,fltry])
        output_filters = np.swapaxes(of,0,1)

        if self.config['VERBOSE']==2:
            print(f'Created square waves at {wavelength} um')

        return output_filters



    def convolve_all_templates(self,idx,norm=True):

        start=time.time()

        z = self.zcat[idx]
        sedx = self.templ*(1+z)
        fl = (1.+z)*(4.*pi*(cosmo.luminosity_distance(z)).to(u.cm).value**2)**-1

        if norm:
            self.dustnorm = (10**10*(const.M_sun/const.m_p)).value


        self.fconv_optical = np.zeros((len(self.sfx),self.optical_grid.shape[1]))
        self.fconv_agn = np.zeros((len(self.sfx),self.agn_grid.shape[1]))
        self.fconv_ir = np.zeros((len(self.sfx),self.ir_grid.shape[1],self.ir_grid.shape[2],self.ir_grid.shape[3]))

        if self.config['EXTRA_BANDS']:
            extrawav = self.wav_extra[idx]
            filt_extra=np.array(self.make_extra_filters(extrawav,5,0.2))
            self.fconv_optical_extra = np.zeros((self.wav_extra.shape[1],self.optical_grid.shape[1]))
            self.fconv_agn_extra = np.zeros((self.wav_extra.shape[1],self.agn_grid.shape[1]))
            self.fconv_ir_extra = np.zeros((self.wav_extra.shape[1],self.ir_grid.shape[1],self.ir_grid.shape[2],self.ir_grid.shape[3]))

        if self.config['FIT_STELLAR']:
            for i,_ in enumerate(self.optical_grid[0,:]):
                for f,_ in enumerate(self.sfx):
                    self.fconv_optical[f,i] = self.convolver(sedx,self.optical_grid[:,i],self.filters[f][0],self.filters[f][1],mode='opt')

                if self.config['EXTRA_BANDS']:
                    for ff,_ in enumerate(extrawav):
                        self.fconv_optical_extra[ff,i]=self.convolver(sedx,self.optical_grid[:,i],filt_extra[ff][0],filt_extra[ff][1],mode='opt')


        if self.config['FIT_AGN']:
            for i,_ in enumerate(self.agn_grid[0,:]):
                for f,_ in enumerate(self.sfx):
                    self.fconv_agn[f,i] = self.convolver(sedx,self.agn_grid[:,i],self.filters[f][0],self.filters[f][1],mode='ir')

                if self.config['EXTRA_BANDS']:
                    for ff,_ in enumerate(extrawav):
                        self.fconv_agn_extra[ff,i]=self.convolver(sedx,self.agn_grid[:,i],filt_extra[ff][0],filt_extra[ff][1],mode='ir')


        if self.config['FIT_DUST']:
            for i,_ in enumerate(self.ir_grid[0,:,0,0]):
                for j,_ in enumerate(self.ir_grid[0,0,:,0]):
                    for k,_ in enumerate(self.ir_grid[0,0,0,:]):
                        sedy = self.ir_grid[:,i,j,k]*self.dustnorm
                        for f,_ in enumerate(self.sfx):
                            self.fconv_ir[f,i,j,k] = self.convolver(sedx,sedy,self.filters[f][0],self.filters[f][1],mode='ir')

                        if self.config['EXTRA_BANDS']:
                            for ff,_ in enumerate(extrawav):
                                self.fconv_ir_extra[ff,i,j,k]=self.convolver(sedx,sedy,filt_extra[ff][0],filt_extra[ff][1],mode='ir')



        if self.config['EXTRA_BANDS']:
            self.fconv_optical=np.concatenate([self.fconv_optical,self.fconv_optical_extra],axis=0)
            self.fconv_agn=np.concatenate([self.fconv_agn,self.fconv_agn_extra],axis=0)
            self.fconv_ir=np.concatenate([self.fconv_ir,self.fconv_ir_extra],axis=0)
        #if self.config['VERBOSE']==2:

        #if norm:
            #self.dustnorm = (10**10*(const.M_sun/const.m_p)).value
            #self.fconv_ir *= self.dustnorm

        if self.config['VERBOSE']==2:
            print(f'Synthethic Photometry Took {time.time()-start} s')

        return None


    def convolver(self,sedx,sedy,filter_x,filter_y,mode='ir'):
        """
        CONVOLVE THE SED WITH THE FILTER OF CHOICE
        MODE

        Parameters
        --------------------
        sedx: sed wavelength
        sedy: sed flux
        filter_x: filter wavelength
        filter_y: filter throughput
        mode: convolution mode, use "opt" when using optical tempalates, "ir" for rest

        """

        integrator=np.trapz

        if mode == 'ir':
            filter_x = filter_x[filter_y>0.0]
            filter_y = filter_y[filter_y>0.0]
            if len(sedx)!=len(filter_y):
                sedy=np.interp(filter_x,sedx,sedy)

        if mode == 'opt':
           sedy = self.interp_conserve(filter_x, sedx,sedy, left=0, right=0)

        temp_int = integrator(filter_y*sedy/filter_x, filter_x)/integrator(filter_y/filter_x, filter_x)
        return temp_int

    def interp_conserve(self,x, xp, fp, left=0., right=0.):

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

    def fit_object(self,idx,scale_sigma=True,parallel=False,vnorm=False):
        """
        FIT A SINGLE GALAXY

        Parameters
        --------------------
        scale_sigma: Rescale the uncertainties by UNCERT_SCALE value in the config, default is True
        """
        fit_start=time.time()


        fnu = np.copy(self.fnu[idx])
        efnu = np.copy(self.efnu[idx])

        efnu[(fnu/efnu)<=3]*=3

        to_fit = fnu>1e-11

        self.usephot[idx,:] = to_fit


        coeffs_obj = np.zeros((total,self.ir_grid.shape[1],self.ir_grid.shape[2],self.ir_grid.shape[3]))
        chi_obj = np.zeros_like(self.ir_grid[0,:,:,:])

        """
        Masking corrupt/missing entries
        """
    
        fnu = fnu[to_fit]
        efnu = efnu[to_fit]
        wav = self.wav[idx,to_fit]

        if scale_sigma:
            efnu[wav<=25] = np.hypot(efnu[wav<=25],self.config['UNCERT_SCALE']*fnu[wav<=25])
            #efnu = np.hypot(efnu,self.config['UNCERT_SCALE']*fnu)
        
        self.nfilt[idx] = len(fnu) #Number of usable bands per object

        self.convolve_all_templates(idx)

        A = np.vstack((self.fconv_optical.T,self.fconv_agn.T,np.zeros_like(self.fconv_agn[:,0])))


        for i,_ in enumerate(self.fconv_ir[0,:,0,0]):
                for j,_ in enumerate(self.fconv_ir[0,0,:,0]):
                    for k,_ in enumerate(self.fconv_ir[0,0,0,:]):
                        B = np.copy(A)
                        B[-1,:] = self.fconv_ir[:,i,j,k]
                        B = B[:,to_fit]
                        B = B/efnu
                        if np.sum(~np.isfinite(B))>0:
                            continue
                        try:
                            coeffs_obj[:,i,j,k],chi_obj[i,j,k] = sp.optimize.nnls(B.T,fnu/efnu)
                        except Exception as e:
                            print(e)
                            coeffs_obj[:,i,j,k], chi_obj[i,j,k] = np.array([0.]*total), -99

        


        chi_obj = chi_obj**2
        best_idx = np.unravel_index(np.argmin(chi_obj),chi_obj.shape)
        A[-1] = self.fconv_ir[:,best_idx[0],best_idx[1],best_idx[2]]

        self.best_ir_idx[idx,:] = best_idx
        self.best_coeffs[idx] = coeffs_obj[:,best_idx[0],best_idx[1],best_idx[2]]
        self.chi2[idx] = chi_obj[best_idx[0],best_idx[1],best_idx[2]]


       
        if not parallel:
            self.compute_physical_params(idx,vnorm=vnorm)
            self.compute_uncertaintites_simple(idx,coeffs_obj,chi_obj,A)
            
            return None

        if parallel:
           
            output = self.mp_output(idx,coeffs_obj,chi_obj,A)

            return idx, output


        
        




    def compute_physical_params(self,idx,vnorm=False):

        norm = self.normalise_coeffs(idx)

        coeffs_rest = (self.best_coeffs[idx,:12].T*norm).T
        coeffs_norm = self.best_coeffs[idx,:12]/np.sum(self.best_coeffs[idx,:12])
        coeffs_rest = np.array(coeffs_rest)

        if vnorm:
            restV,Lv,templf = self.get_rest_flux(idx,154)

            coeffs_vnorm = self.best_coeffs[idx,:12]*templf
            coeffs_vnorm/= coeffs_vnorm.sum()

            Lv_norm = (coeffs_vnorm*self.optical_param['Lv']).sum()
            Lv_norm *= u.solLum
            Mv = (coeffs_vnorm*self.optical_param['mass']).sum()
            Mv *= u.solMass 
            Mv *= 1./Lv_norm

            self.mstar[idx] = (Mv*Lv).value

        else:
            self.mstar[idx] = coeffs_rest.dot(self.optical_param['mass'])


        self.sfr_optical[idx] = coeffs_rest.dot(self.optical_param['sfr'])
        self.av[idx] = coeffs_norm.dot(self.optical_param['Av'])
        self.mdust[idx] = self.dg_ratio*self.best_coeffs[idx,-1]*self.dustnorm*(const.m_p/const.M_sun)/self.lnu_to_fnu[idx]
        self.lir[idx], self.lir_agn[idx] = self.compute_luminosities(self.zcat[idx],self.best_coeffs[idx],self.best_ir_idx[idx])
        self.lir_tot[idx] = self.lir[idx] + self.lir_agn[idx]
        self.sfr[idx] = self.lir[idx]*1e-10
        self.fagn[idx] = self.lir_agn[idx]/self.lir_tot[idx]
        self.deltaGDR[idx] = self.compute_deltaGDR(idx)
        self.mgas[idx] = self.deltaGDR[idx] * self.mdust[idx]


        self.avu[idx] = self.lir[idx]/(125*self.mdust[idx])
        self.umin[idx] = self.templ_umin[np.int_(self.best_ir_idx)[idx,0]]
        self.qpah[idx] = self.q_indices[np.int_(self.best_ir_idx)[idx,1]]
        self.gamma[idx] = self.templ_gamma[np.int_(self.best_ir_idx)[idx,2]]

        return None


    def normalise_coeffs(self,idx):
        '''
        NORMALISE THE OPTICAL COEFFICIENTS
        '''

        z = self.zcat[idx]
        dl=cosmo.luminosity_distance(z).to(u.cm)

        fnu_units = u.erg/u.s/u.cm**2/u.Hz
        mJy_to_cgs = u.mJy.to(u.erg/u.s/u.cm**2/u.Hz)
        fnu_scl = 10**(-0.4*(self.config['ABZP']-23.9))*1*mJy_to_cgs

        template_fnu_units=(1*u.solLum / u.Hz)
        to_physical = fnu_scl*fnu_units*4*np.pi*dl**2/(1+z)
        to_physical /= (1*template_fnu_units).to(u.erg/u.second/u.Hz)

        return to_physical


    def compute_uncertaintites_simple(self,idx,coeffs_obj,chi_obj,A):
        resampled_coeffs = self.resample_coeffs(idx,A)

        norm = self.normalise_coeffs(idx)

        coeffs_rest_resampled = (resampled_coeffs[:,:12].T*norm).T
        coeffs_norm_resampled = (resampled_coeffs[:,:12].T/np.sum(resampled_coeffs[:,:12],axis=1)).T
        coeffs_rest_resampled = np.array(coeffs_rest_resampled)



        mstar_cov=np.tensordot(coeffs_rest_resampled,self.optical_param['mass'],axes=1)
        sfr_optical_cov=np.tensordot(coeffs_rest_resampled,self.optical_param['sfr'],axes=1)
        av_cov=np.tensordot(coeffs_norm_resampled,self.optical_param['sfr'],axes=1)

        #mu = 1.3
        mu = 11.49
        #red = 1/(self.nfilt[idx]-11)

        delta_chi = np.abs(chi_obj-self.chi2[idx])
        coeff_flat = coeffs_obj.reshape(coeffs_obj.shape[0],-1).T[delta_chi.flatten()<mu]
    
        idx_ir = np.unravel_index(np.argwhere(delta_chi<mu),chi_obj.shape)[2]

        lir_cov = np.zeros(len(idx_ir))
        lir_agn_cov = np.zeros(len(idx_ir))

        for i,_ in enumerate(idx_ir):
            lir_cov[i],lir_agn_cov[i] = self.compute_luminosities(self.zcat[idx],coeff_flat[i],idx_ir[i])

        mdust_cov = self.dg_ratio*coeff_flat[:,-1]*self.dustnorm*(const.m_p/const.M_sun)/self.lnu_to_fnu[idx]
        
        self.e_mstar[idx] = utils.to_sigma(mstar_cov)
        self.e_sfr_optical[idx] = utils.to_sigma(sfr_optical_cov)
        self.e_av[idx] = utils.to_sigma(av_cov)
        self.e_lir[idx] = utils.to_sigma(lir_cov)
        self.e_lir_agn[idx] = utils.to_sigma(lir_agn_cov)
        self.e_lir_tot[idx] = utils.to_sigma(lir_cov+lir_agn_cov)
        self.e_fagn[idx] = utils.to_sigma(lir_agn_cov/(lir_agn_cov+lir_cov))
        self.e_mdust[idx] = utils.to_sigma(mdust_cov)
        self.e_mgas[idx] = self.deltaGDR[idx] * utils.to_sigma(mdust_cov)

        self.e_sfr = self.e_lir*10**-10
        self.e_avu = utils.to_sigma(lir_cov/(125*mdust_cov))

        return None


    def resample_coeffs(self,idx,A):

        start = time.time()

        irdx = np.int_(self.best_ir_idx[idx])

        #A = np.vstack((self.fconv_optical.T,self.fconv_agn.T,np.zeros_like(self.fconv_agn[:,0])))
        #A[-1] = self.fconv_ir[:,irdx[0],irdx[1],irdx[2]]

        bestfit=np.copy(A)

        A/=self.efnu[idx]
        A = A[:,np.bool_(self.usephot[idx])]
        A = A.T 

        nonzero_mask = self.best_coeffs[idx]>0
        nonzeroA = A[:,nonzero_mask]

        cov = np.matrix(np.dot(nonzeroA.T, nonzeroA)).I.A

        try:
            covsam = np.random.multivariate_normal(self.best_coeffs[idx,nonzero_mask],cov,5*10**3)
        except Exception as e:
            print(e)
            resampled_coeffs=np.zeros((len(self.best_coeffs[idx,nonzero_mask]),self.best_coeffs.shape[1]))
            return resampled_coeffs

        incl_val = []

        for val,_ in enumerate(covsam[:,0]):
            if all(h>0 for h in covsam[val,:]):
                incl_val.append(val)

        covsam = covsam[incl_val,:]

        resampled_coeffs=np.zeros((len(covsam[:,0]),self.best_coeffs.shape[1]))
        resampled_coeffs[:,self.best_coeffs[idx]>0]=covsam

        #resampled_chi2 = np.zeros((resampled_coeffs.shape[0]))
        

        #for i,_ in enumerate(resampled_coeffs[:,0]):
        #    resampled_chi2[i] = np.sum( (np.sum((resampled_coeffs[i,:]*bestfit.T),axis=1)-self.fnu[idx])**2*(self.efnu[idx]**2)**-1[np.bool_(self.usephot[idx])])



        return resampled_coeffs
        

    

    def compute_deltaGDR(self,idx,use_own_mass=False,allow_own_mass=True):
        '''
        CALCULATE GAS-TO-DUST RATIO 
        
        See Mannucci+10 and Kewley+08 for more details

         Parameters
        --------------------
        use_own_mass: use the stellar mass computed with compute_physical_params, otherwise the catalogue value will be used

        allow_own_mass: if stellar mass in the catalogue is null, switch to using own stellar mass

        '''

        if use_own_mass:
            mstar = self.mstar[idx]
        else:       
            mstar = self.mstar_cat[idx]
        
            if mstar<=0 and allow_own_mass:
                mstar = self.mstar[idx]


        if mstar>0 and self.sfr[idx]>0:
            mu_032=np.log10(mstar)-0.32*np.log10(self.sfr[idx])
            x=mu_032-10.
            if mu_032<10.2:
                met_KD02=8.90+0.47*x
            if mu_032>10.5:
                met_KD02=9.07
            else:
                met_KD02=8.90+0.39*x-0.2*x**2-0.077*x**3+0.064*x**4

            met_PP04=569.4927-192.51820*met_KD02+21.918360*met_KD02**2-0.8278840*met_KD02**3
            deltaGDR=10**(10.54-0.99*met_PP04)
        else:
            #print('Failed to compute deltaGDR')
            deltaGDR=-99

        return deltaGDR



    def compute_luminosities(self,z,coeffs,irdx):
        """
        COMPUTE LIR AND LAGN FROM BEST FIT TEMPLATES

        Parameters
        --------------------
        z: redshift
        coeffs: coefficients
        irdx: [i,j,k] indices of the ir template
        """
        irdx = np.int_(irdx)
        #fl = (1.+z)*(4.*pi*(cosmo.luminosity_distance(z)).to(u.cm).value**2)**-1

        fnu_to_lsol = ((1.0*u.mJy).to(u.erg/u.s/u.Hz/u.cm**2)*4*pi*cosmo.luminosity_distance(z)**2/(1+z)).to(u.Lsun/u.Hz)

        lnu_ir = coeffs[-1]*self.dustnorm*self.ir_grid[:,irdx[0],irdx[1],irdx[2]] * fnu_to_lsol #* fl
        lnu_agn = coeffs[-3:-1]*self.agn_grid * fnu_to_lsol
        lnu_agn = np.sum(lnu_agn,axis=1)

        ir_range = (self.templ>=8) & (self.templ<=1e3)

        nu=(const.c/(self.templ*u.um)).to(u.Hz)
        lir=np.abs(np.trapz(lnu_ir[ir_range],nu[ir_range]).value)
        lagn=np.abs(np.trapz(lnu_agn[ir_range],nu[ir_range]).value)

        return lir,lagn


    def save_sed(self,obj):
        
        if not self.config['SAVE_SED']:
            return None

        sedloc=self.config['PATH']+'seds/'



        #for idx in range(self.fnu.shape[0]):
        for idx in obj:

            irdx = np.int_(self.best_ir_idx[idx])

            sed_x = self.templ  * (1+self.zcat[idx])
            sed_opt = np.sum(self.best_coeffs[idx,:-3]*self.optical_grid,axis = 1)
            sed_agn = np.sum(self.best_coeffs[idx,-3:-1]*self.agn_grid,axis = 1)
            sed_ir = self.best_coeffs[idx,-1]*self.dustnorm*self.ir_grid[:,irdx[0],irdx[1],irdx[2]]

            sed_tot = sed_opt + sed_agn + sed_ir

            table_sed = Table()
            table_sed['lambda'] = sed_x
            table_sed['templ_opt'] = sed_opt
            table_sed['templ_agn'] = sed_agn
            table_sed['templ_ir'] = sed_ir


            table_sed.write(f'{sedloc}/{self.id[idx]}.fits')

        print('Saved all SEDs')


        return None

    def get_obs_flux(self,idx,filter_x,filter_y):

        irdx = np.int_(self.best_ir_idx[idx])

        sed_x = self.templ  * (1+self.zcat[idx])
        sed_opt = np.sum(self.best_coeffs[idx,:-3]*self.optical_grid,axis = 1)
        sed_agn = np.sum(self.best_coeffs[idx,-3:-1]*self.agn_grid,axis = 1)
        sed_ir = self.best_coeffs[idx,-1]*self.dustnorm*self.ir_grid[:,irdx[0],irdx[1],irdx[2]]

        sed_tot = sed_opt + sed_agn + sed_ir

        obs_flux = self.convolver(sed_x,sed_tot,filter_x,filter_y)

        return obs_flux


    def get_rest_flux(self,idx,filter_id):

        filter_x = self.filters_all[filter_id].wave/1e4
        filter_y = self.filters_all[filter_id].throughput

        lspl=self.filters_all[filter_id].name.split()

        for index,s in enumerate(lspl):
            if 'lambda' in s:
                filter_cen = lspl[index+1]*u.Angstrom
                
        

        irdx = np.int_(self.best_ir_idx[idx])

        sed_x = self.templ
        sed_opt = np.sum(self.best_coeffs[idx,:-3]*self.optical_grid,axis = 1)
        sed_agn = np.sum(self.best_coeffs[idx,-3:-1]*self.agn_grid,axis = 1)
        sed_ir = self.best_coeffs[idx,-1]*self.dustnorm*self.ir_grid[:,irdx[0],irdx[1],irdx[2]]

        sed_tot = sed_opt + sed_agn + sed_ir

        lum_dist = (cosmo.luminosity_distance(self.zcat[idx])).to(u.cm)

        rest_flux = self.convolver(sed_x,sed_tot,filter_x,filter_y)

        lnu = (rest_flux*u.mJy).to(u.erg/u.s/u.Hz/u.cm**2)*(4*np.pi*lum_dist**2)/(1+self.zcat[idx])

        L = (lnu*(const.c/filter_cen)).to(u.L_sun) # Luminosity at chosen filter

        templf = []
        for i in range(self.optical_grid.shape[1]):
            x_ = self.templ
            y_ = self.optical_grid[:,i]

            templf.append(self.convolver(x_,y_,filter_x,filter_y))

        return rest_flux,L,templf


    def show_fit(self,idx,xlim=None,ylim=None,components=True,radio=False,detailed=False):
        

        irdx = np.int_(self.best_ir_idx[idx])

        sed_x = self.templ  * (1+self.zcat[idx])
        sed_opt = np.sum(self.best_coeffs[idx,:-3]*self.optical_grid,axis = 1)
        sed_agn = np.sum(self.best_coeffs[idx,-3:-1]*self.agn_grid,axis = 1)
        sed_ir = self.best_coeffs[idx,-1]*self.dustnorm*self.ir_grid[:,irdx[0],irdx[1],irdx[2]]

        sed_tot = sed_opt + sed_agn + sed_ir

        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(1,1,1)

        if components:
            ax.fill_between(sed_x,0,sed_opt,color='royalblue',alpha=0.2,label='Stellar')
            ax.fill_between(sed_x,0,sed_agn,color='g',alpha=0.2,label='AGN')
            ax.fill_between(sed_x,0,sed_ir,color='maroon',alpha=0.2,label='Dust')

        if detailed:
            #To plot optical templates individually
            for i,tempnorm in enumerate(self.best_coeffs[idx,:-3]):
                if tempnorm==0:
                    continue
                ax.plot(sed_x,tempnorm*self.optical_grid[:,i],label=f'{i}')

        if radio:
            sed_radio = self.radio_sed(idx,sed_x,alpha=-0.75)
            sed_tot += sed_radio
        else:
            sed_radio = np.zeros_like(self.templ)

        ax.plot(sed_x,sed_tot,'k',label='Total',lw=3,alpha=0.6,zorder=10)

        points=((self.fnu[idx]/self.efnu[idx])>=3)


        ax.errorbar(self.wav[idx,points],self.fnu[idx,points],yerr=self.efnu[idx,points],color='red',fmt='s',capsize=5,capthick=1,ms=12,markerfacecolor='white',mew=2,barsabove=True)
        ax.scatter(self.wav[idx,~points],(self.fnu[idx,~points]+3*self.efnu[idx,~points]),marker=r'$\downarrow$',s=300,color='red',zorder=11)

        ax.text(0.02, 0.85,'ID '+str(self.id[idx]), color='k',fontsize=20,transform=ax.transAxes)
        ax.text(0.02, 0.75,r'z = {:.2f}'.format(self.zcat[idx]), color='k',fontsize=20,transform=ax.transAxes)

        ax.set_ylabel(r'$f_{\nu}$ [mJy]',fontsize=25)
        ax.set_xlabel(r'$\lambda_{obs}$ $[\mu m]$',fontsize=25)
        ax.legend(fontsize=12)

        ax.set_ylim(10**-4,10**3)
        ax.set_xlim(.5,10**5.7)

        if xlim is not None:
            xlim = xlim
            ax.set_xlim(xlim[0],xlim[1])

        if ylim is not None:
            ylim = ylim
            ax.set_ylim(ylim[0],ylim[1])

        ax.grid(alpha=0.4)
        ax.loglog()

        fig.tight_layout()
        fig.show()

        #if return_output:
        output = {'lambda':sed_x,'templ_optical':sed_opt,'templ_agn':sed_agn,'templ_ir':sed_ir,
            'pivot':self.wav[idx],'fobs':self.fnu[idx],'efobs':self.efnu[idx],'templ_radio':sed_radio}
            

        return fig,output


    def radio_sed(self,idx,lam,alpha=-0.75,use_own_mass=False,allow_own_mass=True):
        '''
        Calculate the radio contribution:
        Two methods available:
        Delvecchio+21, with Lir and Mstar dependence - delv20
        Delhaize+17, only with Lir dependence - delhz17
        '''
        method = self.config['RADIO_METHOD']

        lum_dist = (cosmo.luminosity_distance(self.zcat[idx])).to(u.cm).value

        if self.lir_tot[idx]<0:
            output = 0*lam
            return output


        if method == 'delv20':
            slope = self.radio_slope_delv

            if use_own_mass:
                mstar = self.mstar[idx]
            if not use_own_mass:       
                mstar = self.mstar_cat[idx]

            output=slope(self.zcat[idx],self.lir_tot[idx],mstar,lum_dist,alpha)*lam**(-alpha)

        elif method=='delhz17':
            slope=self.radio_slope_delhz
            output=slope(self.zcat[idx],self.lir_tot[idx],lum_dist,alpha)*lam**(-alpha)

        return output

    def radio_slope_delhz(self,z,LIR,dl,alpha):
        '''
        Compute radio slope as described in Delhaize+17
        '''
        dl=(dl*u.cm).to(u.m)
        q=2.88*(1+z)**(-0.19)
        L14=(10**(np.log10(LIR*(3.839*10**26)/(3.75*10**12))-q))*u.W/u.Hz
        S3=L14*(1+z)**(alpha+1)*(4*np.pi*dl**2)**-1*(3/1.4)**alpha
        S3=S3.to(u.mJy).value
        A=S3*(10**5)**alpha
        return A

    def radio_slope_delv(self,z,LIR,Mstar,dl,alpha):
        '''
        Compute radio slope as described in Delvecchio+21
        '''
        if Mstar<=0:
            return 0.

        dl=(dl*u.cm).to(u.m)
        A=(1+z)
        B=np.log10(Mstar)-10
        q=2.646*A**(-0.023)-B*0.148
        L14=(10**(np.log10(LIR*(3.839*10**26)/(3.75*10**12))-q))*u.W/u.Hz
        S3=L14*(1+z)**(alpha+1)*(4*np.pi*dl**2)**-1*(3/1.4)**alpha
        S3=S3.to(u.mJy).value
        A=S3*(10**5)**alpha
        return A

    def fit_MBB(self,idx,beta=None,kappa=None,type=0):
        """
        FIT A MODIFIED BLACKBODY FUNCTION
        type: 0 - optically thin MBB, 1 - optically thick MBB, 2 - double optically thin MBB
        """

        if beta is None:
            beta = 1.8

        if kappa is None:
            kappa = 0.51*(250/850)**beta


        fnu = np.copy(self.fnu[idx])
        efnu = np.copy(self.efnu[idx])
        restwav = self.wav[idx]/(1+self.zcat[idx])

        efnu[(fnu/efnu)<=3]*=3

        to_fit = (fnu>1e-11) & (restwav>40)

        fnu = fnu[to_fit]
        efnu = efnu[to_fit]





        return None


    def assign_z_bin(self,verbose=False,gridsize=1000,tol=0.01):
        
        '''
        Experimental.
        Boost the fitting efficiency by grouping objects in the catalogue to some common
        redshift. 
        
        Use tol to define the offset tolerance, in units of abs(z_true - zgrid)/(1+z_true)
        values above 0.02 are generally not recommended
        '''


        self.zgrid = np.linspace(np.min(self.zcat),np.max(self.zcat)+0.01,gridsize)

        self.iz = np.zeros(self.NOBJ)
        self.z_uncert = np.zeros(self.NOBJ)
        self.lum_corr = np.zeros(self.NOBJ)

        for i,z in enumerate(self.zcat):
            idz = utils.find_nearest(self.zgrid,z)
            self.iz[i] = idz
            self.z_uncert[i] = np.abs(z-self.zgrid[idz])/(1+z)
            self.lum_corr[i] = (cosmo.luminosity_distance(self.zgrid[idz])**2/(1+self.zgrid[idz]))/(cosmo.luminosity_distance(z)**2/(1+z))

        self.iz = np.int_(self.iz)

        if np.sum(np.abs(self.z_uncert)>tol)!=0:
            print(f'Found a redshift offset in the new grid of more than {tol}, create a finer grid')

        if verbose:
            print(f'Assigned galaxies to a new redshift grid')
            print('Expected performance boost is {:.1%}'.format( (self.NOBJ/len(np.unique(self.iz))-1)))
            print(f'Largest z uncertainty is {np.round(np.max(self.z_uncert),4)}')
            


        return None


    def mp_output(self,idx,coeffs_obj,chi_obj,A):
        
        mp_out = [self.usephot[idx,:],self.nfilt[idx],self.best_ir_idx[idx,:],
                   self.best_coeffs[idx],self.chi2[idx],coeffs_obj,chi_obj,A]

        
        return mp_out

    def mp_restructure(self,mp_out,vnorm=False):

        failed = []
        for i in range(len(mp_out)):
            idx = mp_out[i][0]

            self.usephot[idx,:] = mp_out[i][1][0]
            self.nfilt[idx] = mp_out[i][1][1]
            self.best_ir_idx[idx,:] = mp_out[i][1][2]
            self.best_coeffs[idx] = mp_out[i][1][3]
            self.chi2[idx] = mp_out[i][1][4]
            coeffs_obj = mp_out[i][1][5]
            chi_obj = mp_out[i][1][6]
            A = mp_out[i][1][7]

            try:
                self.compute_physical_params(idx,vnorm=vnorm)

                self.compute_uncertaintites_simple(idx,coeffs_obj,chi_obj,A)
            except:
                failed.append(idx)
                continue
        
        if len(failed)>0:
            print('Failed to get physical params for the following objects: ')
            for idx in failed:
                print(self.id[idx])
        return None

    def fit_catalogue(self, n_proc=-1,n_obj=None,save_results=True,vnorm=False):
        '''
        Main function to fit the full catalogue
        Outputs the ctf class object with all the parameters
         Parameters
        --------------------
        n_proc: number of threads to use, -1 will use all available resources, 1 will fit in serial mode
        '''

        if n_obj is not None:
            obj = range(n_obj)
        else:
            obj = range(self.NOBJ)

        if old_ver:
            import multiprocessing as mp
            global mp_wrapper

        if not old_ver:
            import multiprocess as mp


        def mp_wrapper(i):
            idx,otp = self.fit_object(i,parallel=True)
            return idx,otp


        if n_proc>mp.cpu_count():
                print(f'Requested {n_proc} threads, but only {mp.cpu_count()} are available')


        t0 = time.time()

        #Serial Mode
        
        if (mp.cpu_count()==1) | (n_proc==1):
            print('Running in Serial Mode')
            if HAS_TQDM:
                for idx in tqdm(obj):
                    _res = self.fit_object(idx)
            else:
                for idx in obj:
                    _res = self.fit_object(idx)

            print ('Time Elapsed:',time.time()-t0, 's')
            self.save_results()

            return None


        print('Begin Multithreading')
        if n_proc==-1:
            pool = mp.Pool()
        else:
            pool = mp.Pool(n_proc)

        print(mp.cpu_count(),'threads utilised')

        if HAS_TQDM:
            jobs = [pool.apply_async(mp_wrapper,(i,) ) for i in obj]
            pool.close()
            mp_out = []
            for res in tqdm(jobs):
                mp_out.append(res.get())
            mp_out=np.array(mp_out,dtype=object)

        else:
            mp_out = np.array(pool.map(mp_wrapper,obj))

        print('Finished Fitting, Preparing Output..')

        if vnorm:
            print('Using V-band luminosity to normalise physical parameters')

        self.mp_restructure(mp_out,vnorm=vnorm)

        if save_results:
            self.save_results()

            if self.config['SAVE_FIGURE']:
                self.save_all_figures(obj)

            if self.config['SAVE_SED']:
                self.save_sed(obj)
        print ('Time Elapsed:',time.time()-t0, 's')


        return None


    def save_all_figures(self,obj,xlim=None,ylim=None,components=True,radio=False,detailed=False):

        if self.config['SAVE_FIGURE']:
            figloc=self.config['PATH']+'figures/'
        
        for idx in obj:
            
            irdx = np.int_(self.best_ir_idx[idx])

            sed_x = self.templ  * (1+self.zcat[idx])
            sed_opt = np.sum(self.best_coeffs[idx,:-3]*self.optical_grid,axis = 1)
            sed_agn = np.sum(self.best_coeffs[idx,-3:-1]*self.agn_grid,axis = 1)
            sed_ir = self.best_coeffs[idx,-1]*self.dustnorm*self.ir_grid[:,irdx[0],irdx[1],irdx[2]]

            sed_tot = sed_opt + sed_agn + sed_ir

            fig = plt.figure(figsize=(10,5))
            ax = fig.add_subplot(1,1,1)

            if components:
                ax.fill_between(sed_x,0,sed_opt,color='royalblue',alpha=0.2,label='Stellar')
                ax.fill_between(sed_x,0,sed_agn,color='g',alpha=0.2,label='AGN')
                ax.fill_between(sed_x,0,sed_ir,color='maroon',alpha=0.2,label='Dust')

            if radio:
                sed_radio = self.radio_sed(idx,sed_x,alpha=-0.75)
                sed_tot += sed_radio

            ax.plot(sed_x,sed_tot,'k',label='Total',lw=3,alpha=0.6,zorder=10)

            points=((self.fnu[idx]/self.efnu[idx])>=3)


            ax.errorbar(self.wav[idx,points],self.fnu[idx,points],yerr=self.efnu[idx,points],color='red',fmt='s',capsize=5,capthick=1,ms=12,markerfacecolor='white',mew=2,barsabove=True)
            ax.scatter(self.wav[idx,~points],(self.fnu[idx,~points]+3*self.efnu[idx,~points]),marker=r'$\downarrow$',s=300,color='red',zorder=11)

            ax.text(0.02, 0.85,'ID '+str(self.id[idx]), color='k',fontsize=20,transform=ax.transAxes)
            ax.text(0.02, 0.75,r'z = {:.2f}'.format(self.zcat[idx]), color='k',fontsize=20,transform=ax.transAxes)

            ax.set_ylabel(r'$f_{\nu}$ [mJy]',fontsize=25)
            ax.set_xlabel(r'$\lambda_{obs}$ $[\mu m]$',fontsize=25)
            ax.legend(fontsize=12)

            ax.set_ylim(10**-4,10**3)
            ax.set_xlim(.5,10**5.7)

            if xlim is not None:
                xlim = xlim
                ax.set_xlim(xlim[0],xlim[1])

            if ylim is not None:
                ylim = ylim
                ax.set_ylim(ylim[0],ylim[1])

            ax.grid(alpha=0.4)
            ax.loglog()

            fig.savefig(f'{figloc}/{self.id[idx]}.pdf')
            plt.close()

        
                
        print('Saved all figures')

        return None

        

    def save_results(self):

        #Ignore this list
        outputnames=['id','z','mstar_input',
                        'lir_tot','e_lir_tot','lir','e_lir','lir_agn','e_lir_agn', #Luminosity
                        'f_agn','e_f_agn', #AGN
                        'mdust','e_mdust','mgas','e_mgas', 'delta_gdr', #IR masses
                        'avu','e_avu', #IR luminosity/mass
                        'mstar','e_mstar','av','e_av','sfr_opt','e_sfr_opt', #Optical masses and other properties
                        'Umin','qpah','gamma', #Best fit IR params
                        'chi2','nfilt','lastdet', #Misc 
        ]

        self.tab = Table()

        #External
        self.tab['id'] = self.id
        self.tab['ra'] = self.ra
        self.tab['dec'] = self.dec
        self.tab['z'] = self.zcat
        self.tab['mstar_input'] = self.mstar_cat

        #Luminosity
        self.tab['lir'] = self.lir 
        self.tab['e_lir'] = self.e_lir
        self.tab['lir_agn'] = self.lir_agn
        self.tab['e_lir_agn'] = self.e_lir_agn
        self.tab['lir_tot'] = self.lir_tot
        self.tab['e_lir_tot'] = self.e_lir_tot

        self.tab['sfr'] = self.sfr
        self.tab['e_sfr'] = self.e_sfr

        #AGN
        self.tab['f_agn'] = self.fagn
        self.tab['e_f_agn'] = self.e_fagn

        #IR Masses
        self.tab['mdust'] = self.mdust
        self.tab['e_mdust'] = self.e_mdust
        self.tab['mgas'] = self.mgas
        self.tab['e_mgas'] = self.e_mgas
        self.tab['delta_gdr'] = self.deltaGDR

        #Optical Masses
        self.tab['mstar'] = self.mstar
        self.tab['e_mstar'] = self.e_mstar
        self.tab['av'] = self.av
        self.tab['e_av'] = self.e_av
        self.tab['sfr_opt'] = self.sfr_optical
        self.tab['e_sfr_opt'] = self.e_sfr_optical

        #DL07 specific
        self.tab['avu'] = self.avu
        self.tab['e_avu'] = self.e_avu
        self.tab['umin'] = self.umin
        self.tab['qpah_idx'] = self.qpah
        self.tab['gamma'] = np.int_(self.gamma)

        self.tab['chi2'] = self.chi2
        self.tab['nfilt'] = np.int_(self.nfilt)
        self.tab['lastdet'] = self.lastdet

        if self.config['SAVE_TABLE']:
            table_out = self.config['PATH']+self.config['OUTPUT_NAME']+'.fits'
            if exists(table_out):
                print(f'Output table {table_out} already exists')
                print('New output table was not saved, delete the old version and run "self.save_results()" again ')    
            else:    
                self.tab.write(table_out)
        
        return None



    def set_zgrid(self):
        'Setup the redshift grid for photoz calculation'

        try:
            zrange = [self.config['Z_MIN'],self.config['Z_MAX']]
            zstep = self.config['Z_STEP']
        except:
            print('Z_MIN, Z_MAX or Z_STEP parameters are not defined in config.')
            print('Using default values (0.01,12,0.005)')
            zrange = [0.01,12]
            zstep = 0.005

        assert zrange[0]>0, 'Z_MIN IS ZERO OR BELOW ZERO'

        self.zgrid = np.arange(*zrange,zstep)

        return None

    def set_ztemps(self):

        self.set_zgrid()

        fit_range = self.sfx<10.
        filts = []
        for i,wl in enumerate(self.sfx):
            if wl<10.:
                filts.append(self.filters[i])


        self.fconv_optical_zgrid = np.zeros((len(self.sfx[fit_range]),len(self.zgrid),self.optical_grid.shape[1]))

        #MULTIPROCESS THIS

        st = time.time()
        for j,z_i in enumerate(self.zgrid):
            sedx = self.templ*(1+z_i)
            for i,_ in enumerate(self.optical_grid[0,:]):
                for f,_ in enumerate(self.sfx[fit_range]):
                    self.fconv_optical_zgrid[f,j,i] = self.convolver(sedx,self.optical_grid[:,i],filts[f][0],filts[f][1],mode='ir')
        print(f'Redshift Templates Ready. Time Elapsed:{time.time()-st}s')

        return None



    def fit_photoz(self,idx,scale_sigma=True,parallel=False):
        """
        FIT PHOTOZ FOR A SINGLE GALAXY

        Parameters
        --------------------
        scale_sigma: Rescale the uncertainties by UNCERT_SCALE value in the config, default is True
        """
        fit_start=time.time()


        fnu = np.copy(self.fnu[idx])
        efnu = np.copy(self.efnu[idx])

        efnu[(fnu/efnu)<=3]*=3

        to_fit = fnu>1e-11

        #self.usephot[idx,:] = to_fit


        coeffs_obj = np.zeros((total,self.ir_grid.shape[1],self.ir_grid.shape[2],self.ir_grid.shape[3]))
        chi_obj = np.zeros_like(self.ir_grid[0,:,:,:])

        """
        Masking corrupt/missing entries
        """
    
        fnu = fnu[to_fit]
        efnu = efnu[to_fit]
        wav = self.wav[idx,to_fit]

        if scale_sigma:
            efnu[wav<=25] = np.hypot(efnu[wav<=25],self.config['UNCERT_SCALE']*fnu[wav<=25])
            #efnu = np.hypot(efnu,self.config['UNCERT_SCALE']*fnu)
        
        
        self.nfilt_zphot[idx] = len(fnu) #Number of usable bands per object

        self.convolve_all_templates(idx)

        A = np.vstack((self.fconv_optical.T,self.fconv_agn.T,np.zeros_like(self.fconv_agn[:,0])))


        for i,_ in enumerate(self.fconv_ir[0,:,0,0]):
                for j,_ in enumerate(self.fconv_ir[0,0,:,0]):
                    for k,_ in enumerate(self.fconv_ir[0,0,0,:]):
                        B = np.copy(A)
                        B[-1,:] = self.fconv_ir[:,i,j,k]
                        B = B[:,to_fit]
                        B = B/efnu
                        if np.sum(~np.isfinite(B))>0:
                            continue
                        try:
                            coeffs_obj[:,i,j,k],chi_obj[i,j,k] = sp.optimize.nnls(B.T,fnu/efnu)
                        except Exception as e:
                            print(e)
                            coeffs_obj[:,i,j,k], chi_obj[i,j,k] = np.array([0.]*total), -99

        


        chi_obj = chi_obj**2
        best_idx = np.unravel_index(np.argmin(chi_obj),chi_obj.shape)
        A[-1] = self.fconv_ir[:,best_idx[0],best_idx[1],best_idx[2]]

        self.best_ir_idx[idx,:] = best_idx
        self.best_coeffs[idx] = coeffs_obj[:,best_idx[0],best_idx[1],best_idx[2]]
        self.chi2[idx] = chi_obj[best_idx[0],best_idx[1],best_idx[2]]


       
        if not parallel:
            self.compute_physical_params(idx,vnorm=vnorm)
            self.compute_uncertaintites_simple(idx,coeffs_obj,chi_obj,A)
            
            return None

        if parallel:
            output = self.mp_output(idx,coeffs_obj,chi_obj,A)

            return idx, output

    def apply_zeropoint(self):
        zp = np.zeros_like(self.fnu[0,:])
        try:
            lines = open(self.zp_file).readlines()
            for line in lines:
                if not line.startswith('F'):
                    continue
            
                fnum = int(line.strip().split()[0][1:])
                if str(fnum) in self.f_numbers:
                    ix = np.int_(self.f_numbers) == fnum
                    zp[ix] = float(line.split()[1])
            print(f'Added zeropoint corrections from ' + self.zp_file)

            return zp
        except:
            print('Zeropoint file not defined. Skipping...')

            return None



