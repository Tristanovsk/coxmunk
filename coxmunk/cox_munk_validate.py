import coxmunk
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmocean as cm

import coxmunk.coxmunk as coxmunk
from coxmunk.visu import plot as uplot

wind = 14.5
wind_azi = 75
stats =['cm_iso','cm_dir','bh2006']
stat=stats[2]
shadow=True
refrac_index=1.334

sza = 39
vza = np.linspace(0,80,81)
azi = np.linspace(0,360,181)
Nvza,Nazi =len(vza),len(azi)

data=np.zeros((Nazi,Nvza,3))
for i in range(Nvza):
    for j in range(Nazi):
        data[j,i,:]=coxmunk.sunglint(
            sza,vza[i],azi[j],m=refrac_index).sunglint(
            wind, wind_azi, stats=stat, shadow=shadow)

fig, axs = plt.subplots(nrows=2,ncols=2, subplot_kw=dict(projection='polar'),figsize=(15, 13))
axs=axs.ravel()
axs[0].scatter([0],[sza],marker='*',facecolor ='orange',alpha=0.6,s=1000)
uplot().label_polplot(axs[0],yticks=[20., 40., 60.,80.], ylabels=['$20^{\circ}$', '$40^{\circ}$', '$60^{\circ}$',''])

axs[0].arrow(wind_azi*np.pi/180, 0, 0,np.max(vza)*max(3,min(16,wind))/22, alpha = 0.5, width = 0.05,head_width=0.25,head_length=15,
                 edgecolor = 'black', facecolor = 'green', lw = 1.2)
axs[0].set_title('Sun and wind directions',  pad=30)

for i,title in  enumerate(('I','Q','U')):
    if title=='I':
        cmap = plt.cm.gist_stern_r
    else:
        cmap = cm.tools.crop_by_percent(cm.cm.balance, 20, which='both', N=None)
    uplot().add_polplot(axs[i+1], vza, azi, data[...,i].T, title=title, cmap=cmap)
# I=data[...,0].T
# Q=data[...,1].T
# U=data[...,2].T
# DOP = np.sqrt(Q**2+U**2)#/I
# uplot().add_polplot(axs[3], vza, azi, DOP, title='DOP', cmap=cmap)
plt.suptitle(r'Wind speed: {:.1f} m/s; direction: {:.1f} deg.'.format(wind,wind_azi))
plt.tight_layout(rect=[0.0, 0.0, 0.99, 0.95])

# plot Pdist
def Pdist_(ws,eta,xi,stats='bh2006'):
    if stats == 'cm_dir':  # historical values from directional COX MUNK
        s_cr2 = 0.003 + 1.92e-3 * ws
        s_up2 = 3.16e-3 * ws
        s_cr = np.sqrt(s_cr2)
        s_up = np.sqrt(s_up2)
        c21 = 0.01 - 8.6e-3 * ws
        c03 = 0.04 - 33.e-3 * ws
        c40 = 0.40
        c22 = 0.12
        c04 = 0.23
    elif stats == 'bh2006':  # from Breon Henriot 2006 JGR
        s_cr2 = 3e-3 + 1.85e-3 * ws
        s_up2 = 1e-3 + 3.16e-3 * ws
        s_cr = np.sqrt(s_cr2)
        s_up = np.sqrt(s_up2)
        c21 = -9e-4 * ws ** 2
        c03 = -0.45 / (1 + np.exp(7. - ws))
        c40 = 0.30
        c22 = 0.12
        c04 = 0.4
    return np.exp(-5e-1 * (xi ** 2 + eta ** 2)) / (2. * np.pi * s_cr * s_up) * \
                    (1. -
                     c21 * (xi ** 2 - 1.) * eta / 2. -
                     c03 * (eta ** 3 - 3. * eta) / 6. +
                     c40 * (xi ** 4 - 6. * eta ** 2 + 3.) / 24. +
                     c04 * (eta ** 4 - 6. * eta ** 2 + 3.) / 24. +
                     c22 * (xi ** 2 - 1.) * (eta ** 2 - 1.) / 4.)

plt.figure()
ws=10
xi=np.linspace(0,10,100)
eta=1e-2
plt.plot(xi,Pdist_(ws,eta,xi))

