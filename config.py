'INPUT PARAMETERS'
#===============================================================================
cat_name='/Users/vasily/Documents/Tables/sdc+laigle_corrected_mw.fits'

bands_file='/Users/vasily/Documents/PhD/Projects/release/configs/vk/cosmos2015.bands'
param_file='/Users/vasily/Documents/PhD/Projects/release/configs/vk/cosmos2015.param'

extra_bands_file='/Users/vasily/Documents/PhD/Projects/release/configs/vk/cosmos2020_alma.bands_extra'

#===============================================================================

'OUTPUT PARAMETERS'
#===============================================================================
PATH = '/Users/vasily/Documents/PhD/Projects/release/output/chi2_test/'

output_name='chi2_test' #Name of output table

figloc=f'{PATH}figures/' #Location of figures

sedloc=f'{PATH}seds/' #Location of SEDs

covarloc=f'{PATH}covar_tables/' #Location of covariance tables
#===============================================================================

'GENERAL SETTINGS'
#===============================================================================
flux_unit='mJy'
multithread=0
verbose=1 # Allowed 0,1,2
extra_bands=0

use_cold_dl=1

diagplot=1
radio=1
radio_method='delv20'

save_fig=0
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
igm_switch=1
use_own_stellar_mass=1
ABZP=23.9
#===============================================================================
'TEMPLATE PARAMETERS'
#===============================================================================
'SWITCH ON/OFF different components'
dust_switch=1
agn_switch=1
stellar_switch=1
#===============================================================================
