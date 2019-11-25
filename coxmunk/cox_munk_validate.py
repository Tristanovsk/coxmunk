import coxmunk
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lut as l

sza = 45
vza = np.linspace(0,80,73)
azi = np.linspace(0,360,73)
Nvza,Nazi =len(vza),len(azi)

data=np.zeros((Nazi,Nvza))
for i in range(Nvza):
    for j in range(Nazi):
        data[j,i]=coxmunk.sunglint(sza,vza[i],azi[j],m=1.334).sunglint(1.6,41,stats="bh2006")

l.lut().plot_lut(vza,azi,data)

