
'INPUT PARAMETERS'
#===============================================================================
#cat_name='/Users/vasily/Documents/Tables/COSMOS_2020_FIR.fits'
cat_name='/Users/vasily/Documents/PhD/Projects/release/ctf/example/COSMOS_SDC2_example.fits'

bands_file='/Users/vasily/Documents/PhD/Projects/release/ctf/example/example.bands'
param_file='/Users/vasily/Documents/PhD/Projects/release/ctf/example/example.param'

extra_bands_file=''

#===============================================================================

'OUTPUT PARAMETERS'
#===============================================================================
PATH = '/Users/vasily/Documents/PhD/Projects/release/ctf/example/output/'

output_name='example_output' #Name of output table

figloc=f'{PATH}figures/' #Location of figures

sedloc=f'{PATH}seds/' #Location of SEDs

covarloc=f'{PATH}covar_tables/' #Location of covariance tables
#===============================================================================

'GENERAL SETTINGS'
#===============================================================================
multithread=True
verbose=2 # Allowed 0,1,2
extra_bands=False

use_cold_dl=True

diagplot=False
radio=True

save_fig=False
save_table=True
save_sed=True
save_covar=False

impose_cut=False
#===============================================================================

'TEMPLATE PARAMETERS'
#===============================================================================
'SWITCH ON/OFF different components'
dust_switch=1
agn_switch=1
stellar_switch=1
#===============================================================================








'''ADVANCED SETTING BELOW'''
'''ONLY CHANGE IF YOU KNOW WHAT YOU ARE DOING'''
#===============================================================================
'''Modules'''

from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy import constants as const
import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
import sys


print('------------')
print('INITIALISED')
print('------------')



cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
c=const.c
Lsol=(const.L_sun.cgs).value
Msol=(const.M_sun).value
pi=np.pi

#===============================================================================
'''
LOAD IN DL07 Templates
'''
model1 = fits.open('templates/DL07/model1.fits', memmap=True)
model2 = fits.open('templates/DL07/model2.fits', memmap=True)
lam = fits.open('templates/DL07/lambda.fits', memmap=True)


model1_full = Table.read('templates/DL07/model1_full.fits', memmap=True)
model2_full = Table.read('templates/DL07/model2_full.fits', memmap=True)



if use_cold_dl:
    diffuse=np.array(model1_full['MODEL1'][0]) #DL diffuse part as 3d array
    pdr=np.array(model2_full['MODEL1'][0])  #DL pdr part as 3d array
    dat1=diffuse
    dat2=pdr
    Umin=np.array([0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,1.000,1.200,
    1.700,2.000,3.000,4.000,5.000,7.000,8.000,10.00,12.00,15.00,20.00,25.00,30.00,
    35.00,40.00,50.00])
    q_indices=np.array([0,1,2,3,4,5,6,7,8,9,10])
    U_cut=0
    Umin=Umin[:]
    dat1=dat1[:,:,q_indices]
    dat2=dat2[:,:,q_indices]
else:
    dat1=model1[0].data
    dat2=model2[0].data
    Umin=np.array([0.800,1.200,2.000,3.000,4.000,5.000,7.000,8.000,10.00,12.00,
    15.00,20.00,25.00,30.00,35.00,40.00,50.00])

DG_ratio=0.01 #Template specific dust-to-gas ratio, see DL07 for more details.


dat3=lam[0].data
#===============================================================================

DATA=Table.read(f'{cat_name}', memmap=True)

#DATA=DATA[DATA['ID_FIR']==641953]


print(f'Read a catalogue with {len(DATA)} objects')


if save_table:
    table_out=(f'{PATH}{output_name}.fits')
#===============================================================================

'''DL07 CONFIGURATION PANEL'''

'Range of the gamma parameter, change at your discretion'
g1=np.array([0,0.001,0.0025,0.005,0.0075])
g2=np.arange(0.01,0.1,0.01)
g3=np.array([0.2,0.35,0.5])
g=np.concatenate((g1,g2,g3))

#===============================================================================
bnds=np.loadtxt(bands_file,dtype='str')
prms=np.loadtxt(param_file,dtype='str')

if verbose>0:
    print(f'Read band file: {bands_file}')
    print(f'Read param file: {param_file}')

filtname=bnds[:,0].tolist()
band_names=bnds[:,1].tolist()
err_band_names=bnds[:,2].tolist()

#===============================================================================
from filter_importer import *

filt_id=np.int_(filtname).tolist()

FILTERS=[]
sfx=[]
names=[]
for id in filt_id:
    fltr=np.array([filters[id-1].wave*10**-4,filters[id-1].throughput])
    FILTERS.append(fltr)
    lspl=filters[id-1].name.split()
    names.append(lspl[0])
    for index,s in enumerate(lspl):
        if 'lambda' in s:
            sfx.append(lspl[index+1])


sfx=10**-4*np.float_(sfx)

filter_diagnostic=False
if filter_diagnostic:
    import math
    xlen=6
    ylen=math.ceil(len(FILTERS)/6)
    fig=figure(figsize=(3*xlen,3*ylen))

    for f,_ in enumerate(FILTERS[:]):
        print()
        ax = fig.add_subplot(ylen,xlen,f+1)
        ax.plot(FILTERS[f][0],FILTERS[f][1])
        ax.axvline(sfx[f],color='r')
        ax.text(0.2,0.2,np.round(sfx[f],2),transform=ax.transAxes)

    show()
    sys.exit()
#===============================================================================
'''Getting Parameters'''

param_names=prms[:,1].tolist()

if param_names[2]=='None':
    DATA['Mstar']=-99
    param_names[2]='Mstar'
#===============================================================================
assert len(filtname)==len(band_names)==len(err_band_names),'ERROR: Filter range does not match available photometry, check band names'

if verbose>0:
    if verbose>1:
        print(f'Detected the following filters...')
        print('--------------------------------')
        for name in names:
            print(name)
        print('--------------------------------')
    print(f'Detected the following parameters {param_names}')
    print('Input Bands: ok')
    print('Input Filters: ok')



'If your dataset contains some exotic bands, i.e. continuum from a line input the information here'
'You also have to specify the observed wavelength columns'

if extra_bands:
    xtra_bnds=np.loadtxt(extra_bands_file,dtype='str')
    bands_extra=xtra_bnds[:,0].tolist()
    err_bands_extra=xtra_bnds[:,1].tolist()
    wavelength_extra=xtra_bnds[:,2].tolist()
    if verbose==1:
        print(f'Added extra bands from {extra_bands_file}')
#===============================================================================

steltemp=12 #Default value
agntemp=2 #Default value
irtemp=1
total=steltemp+agntemp+irtemp

detection_threshold=3 #In sigma
ir_detections=3 #
ir_cutoff=20.
#===============================================================================

if impose_cut:
    zphot=DATA['zphot']
    ztype=np.ones(len(zphot))
    zbest=DATA['zspec']
    zspec=DATA['zspec']
    ztype[zbest<0]=0
    zbest[zbest<0]=zphot[zbest<0]
    DATA['zphot']=zbest
else:
    ztype=np.array([1]*len(DATA))


#===============================================================================

G=np.array((np.array(DATA[band_names])).tolist())
E=np.array((np.array(DATA[err_band_names])).tolist())
W=np.tile(sfx,(len(G[:,0]),1))
if extra_bands:
    G_ex=np.array((np.array(DATA[bands_extra])).tolist())
    E_ex=np.array((np.array(DATA[err_bands_extra])).tolist())
    W_ex=np.array((np.array(DATA[wavelength_extra])).tolist())
    G=np.concatenate((G,G_ex),axis=1)
    E=np.concatenate((E,E_ex),axis=1)
    W=np.concatenate((W,W_ex),axis=1)


P=np.array((np.array(DATA[param_names])).tolist())

#===============================================================================

if impose_cut:
    SN=G[:,:]/E[:,:]
    mask=[]
    ir_det=[]
    for i,obj in enumerate(SN[:,0]):
        a=SN[i,:]
        b=W[i,:]/(1+zbest[i])
        a=a[b>=ir_cutoff]
        det3=len(a[a>=detection_threshold])
        if det3>=ir_detections and P[i,1]>0:
            #ir_det.append(det3)
            mask.append(i)
    G=G[mask]
    E=E[mask]
    P=P[mask]
    W=W[mask]
    print(f'Specified {ir_detections} detections beyond {ir_cutoff} microns rest frame')

#===============================================================================

print('Data loaded')
print('Objects to fit:',len(G[:,0]))

#===============================================================================

if save_fig:
    try:
        os.mkdir(figloc)
    except:
        pass

if save_sed:
    try:
        os.mkdir(sedloc)
    except:
        pass

if save_covar:
    try:
        os.mkdir(covarloc)
    except:
        pass


if save_table:
    outputnames=['ID','LIR_total','eLIR_total','MD',
    'eMD','z','chi2','f_agn','efagn','lastdet','ztype','MG','eMG',
    'deltaGDR','attempts','Mstar','fgas','fgas_FMR','Lir_med',
    'eLir68','Mdust_med','eMdust68','Umin','qpah','gamma','U','sU',
    'Lagn','eLagn','Lir_draine','eLir_draine']

    dtype=[float]*len(outputnames)
    dtype[0]=int

    table_0=Table(np.zeros(len(outputnames)), names=outputnames,dtype=dtype) #Make empty table
    table_0.write(table_out,overwrite=True)

#===============================================================================
print('End Config')
