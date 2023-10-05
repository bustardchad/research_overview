import numpy as np
import matplotlib.pyplot as plt
# Athena++ history data
# [1]=time     [2]=dt       [3]=mass     [4]=1-mom    [5]=2-mom    [6]=3-mom    [7]=1-KE     [8]=2-KE     [9]=3-KE     [10]=tot-E   [11]=1-ME    [12]=2-ME    [13]=3-ME  
time, ke1, ke2, ke3, tote, me1, me2, me3  = np.loadtxt('cr.hst',skiprows = 2, usecols = (0,6,7,8,9,10,11,12),unpack=True)

edenstocgs = 6.54e-11
#cellVol = (250*3.0856e18)**3.0
cellVol = 0.25**3.0
ketot = []
metot = []
thermaletot = []
for j in range(len(time)):
  keval = np.sqrt(ke1[j]**2 + ke2[j]**2 + ke3[j]**2)*edenstocgs/cellVol
  meval = np.sqrt(me1[j]**2 + me2[j]**2 + me3[j]**2)*edenstocgs/cellVol
  thermalval = (tote[j]*edenstocgs/cellVol) - keval - meval
  ketot.append(keval)
  metot.append(meval)
  thermaletot.append(thermalval)
  
plt.plot(time,thermaletot,'k-',label = r'E$_{th}$',linewidth=3)
plt.plot(time,ketot,'r-',label = r'E$_{KE}$',linewidth=3)
plt.plot(time,metot,'g-',label = r'E$_{B}$',linewidth=3)
plt.xlabel('Time (Myrs)',fontsize=18)
plt.ylabel('Energy Density',fontsize=18)
plt.legend()
plt.savefig('EnergyDensity.pdf')
plt.close()


