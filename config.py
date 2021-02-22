'INPUT PARAMETERS'
#===============================================================================
cat_name='/Users/vasily/Documents/Tables/weirdguys_cls.fits'

bands_file='/Users/vasily/Documents/PhD/Projects/release/configs/vk/cosmos2020_jw.bands'
param_file='/Users/vasily/Documents/PhD/Projects/release/configs/vk/cosmos2020_jw.param'

extra_bands_file=''

#===============================================================================

'OUTPUT PARAMETERS'
#===============================================================================
PATH = '/Users/vasily/Documents/PhD/Projects/release/output/jw_test/'

output_name='jw_220221' #Name of output table

figloc=f'{PATH}figures/' #Location of figures

sedloc=f'{PATH}seds/' #Location of SEDs

covarloc=f'{PATH}covar_tables/' #Location of covariance tables
#===============================================================================

'GENERAL SETTINGS'
#===============================================================================
flux_unit=''
multithread=1
verbose=1 # Allowed 0,1,2
extra_bands=0

use_cold_dl=1

diagplot=1
radio=1

save_fig=1
save_table=1
save_sed=1
save_covar=0

impose_cut=0
#===============================================================================

'TEMPLATE PARAMETERS'
#===============================================================================
'SWITCH ON/OFF different components'
dust_switch=1
agn_switch=1
stellar_switch=1
#===============================================================================
