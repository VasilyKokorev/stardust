import idlsave
import numpy as np
import pandas as pd
from astropy.table import Table
#ss=idlsave.read('survey.save',python_dict=True)
#s=idlsave.read('survey.save',python_dict=False)

#output=np.zeros((len(s.id),52))
#names=[]
#for i in ss:
#   names.append(i)


s = idlsave.read('model1_full.save')
df=pd.DataFrame(s)
as_tab=Table.from_pandas(df)
as_tab.write('model1_full.fits')

'''
for i,_ in enumerate(output[:,0]):
    output[i,:]=[s.id[i],s.ra[i],s.de[i],s.zsp[i],s.f36[i],s.ef36[i],s.f45[i],s.ef45[i]
                 s.f58[i],s.ef58[i],s.f8[i],s.ef8[i],s.f24[i],s.ef24[i],s.f100[i],s.ef100[i]
                 s.f160[i],s.ef160[i],s.f250[i],s.ef250[i],s.f350[i],s.ef350[i],s.f500[i]
                 s.ef500[i],s.f850[i],s.ef850[i],s.f1100[i],s.ef1100[i],s.f1200[i],s.ef1200[i]
                 s.radio[i],s.eradio[i],s.sn_ir[i],s.mstar[i],s.zz[i],s.f3g[i],s.ef3g[i],s.f14g[i],
                 s.ef14g[i],s.f_a13[i],s.ef_a13[i],s.l_a13[i],s.f_a11[i],s.ef_a11[i],s.l_a11[i],
                 s.f_a30[i],s.ef_a30[i],s.l_a30[i],s.sed_id_shuowen[i],s.f_a9[i],s.ef_a9[i],s.l_a9[i]]
'''

