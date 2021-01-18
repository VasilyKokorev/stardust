'INPUT PARAMETERS'
#===============================================================================
cat_name='/example/COSMOS_SDC2_example.fits'

bands_file='/example/example.bands'
param_file='/example/example.param'

extra_bands_file=''

#===============================================================================

'OUTPUT PARAMETERS'
#===============================================================================
PATH = '/example/output/'

output_name='example' #Name of output table

figloc=f'{PATH}figures/' #Location of figures

sedloc=f'{PATH}seds/' #Location of SEDs

covarloc=f'{PATH}covar_tables/' #Location of covariance tables
#===============================================================================

'GENERAL SETTINGS'
#===============================================================================
flux_unit=''
multithread=True
verbose=1 # Allowed 0,1,2
extra_bands=False

use_cold_dl=True

diagplot=True
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
