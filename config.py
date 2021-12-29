'INPUT PARAMETERS'
#===============================================================================
cat_name='example/stellar+ir/c2015_sdc2_example.fits'

bands_file='example/stellar+ir/example.bands'
param_file='example/stellar+ir/example.param'



extra_bands_file=''

#===============================================================================

'OUTPUT PARAMETERS'
#===============================================================================
PATH = 'example/stellar+ir/output/'

output_name='example' #Name of output table

figloc=f'{PATH}figures/' #Location of figures

sedloc=f'{PATH}seds/' #Location of SEDs

covarloc=f'{PATH}covar_tables/' #Location of covariance tables
#===============================================================================

'GENERAL SETTINGS'
#===============================================================================
flux_unit='mJy'
multithread=1
verbose=1 # Allowed 0,1,2
extra_bands=0

use_cold_dl=1

diagplot=1
radio=1
radio_method='delv20'

save_fig=1
save_table=1
save_sed=1
save_covar=0

impose_cut=0
impose_detection_cut=0
#===============================================================================

'ADVANCED SETTINGS'
#===============================================================================
uncert_scale=0.05
qso=0
igm_switch=0 #only works if EAZY is installed
use_own_stellar_mass=0
ABZP=23.9
#===============================================================================
'TEMPLATE PARAMETERS'
#===============================================================================
'SWITCH ON/OFF different components'
dust_switch=1
agn_switch=1
stellar_switch=1
#===============================================================================
