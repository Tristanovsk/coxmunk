
import numpy as np
import xarray as xr
import pandas as pd
import cmocean as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.font_manager import FontProperties

import coxmunk.coxmunk as coxmunk

slope=None
shadow=False
stats='cm_dir'
#stats='cm_iso'
title=""
if stats=='cm_dir':
    title = 'Cox-Munk anisotropic'
elif stats=='cm_iso':
    title = 'Cox-Munk isotropic'

#data = np.zeros((Nazi, Nvza, N))
azi=180
vzas = np.linspace(0, 80, 81)
azis = np.array([0,45,90,135,180])
winds = np.array([0.5,2,6,10,14])
wind_azi=0
szas = np.arange(10,90,10)

xr.Dataset()
arr = []
for wind in winds:
    for sza in szas:
        for azi in azis:
            for vza in vzas:
                I,Q,U,V = coxmunk.sunglint(sza, vza, azi, m=1.334).sunglint(
                    wind, wind_azi, stats=stats, shadow=shadow, slope=slope)
                arr.append([wind,wind_azi,sza,vza,azi, I,Q,U,V])
nparr = np.array(arr)

df = pd.DataFrame(nparr,columns=['wind','wind_azi','sza','vza','azi','I','Q','U','V']).set_index(['wind','wind_azi','sza','vza','azi'])

xarr = df.to_xarray()

xarr_ =xarr.sel(azi=180).squeeze()

xarr_.I.plot(hue='wind',col='sza',col_wrap=4,sharey=False)
plt.suptitle(title)
plt.savefig('test/transect/'+title+'azi{:.1f}'.format(azi)+'.png',dpi=300)