import yt
import trident
from trident import LightRay
import aplpy
import numpy as np
#import matplotlib as plt
from yt.units.yt_array import YTQuantity
from yt import YTArray
import h5py
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import wcs
from astropy.wcs import WCS
import csv
import pandas as pd
import warnings









inputArr = np.genfromtxt('aligned_LMC_SFR_RP_rotation_ChadVersion.csv',delimiter=',',invalid_raise=False)
SFR = inputArr[:,0]
SFR = SFR[1:len(SFR)]
ra = inputArr[:,1]
ra = ra[1:len(ra)]
dec = inputArr[:,2]
dec = dec[1:len(dec)]
RP = inputArr[:,3]
RP = RP[1:len(RP)]

plt.plot(RP, SFR,'bo')
plt.title("RP vs SFR (rotation included)")
plt.xlabel("Ram Pressure",fontsize=16)
plt.ylabel("SFR",fontsize=16)
plt.tight_layout()
plt.savefig("SFR_RP_wholeDisk.pdf")
plt.close()


# cut into quadrants:
Dor30_SFR = SFR[dec < -68]
print(Dor30_SFR)
Dor30_ra = ra[dec < -68]
Dor30_SFR = Dor30_SFR[Dor30_ra > 80]

Dor30_RP = RP[dec < -68]
Dor30_ra = ra[dec < -68]
Dor30_RP = Dor30_RP[Dor30_ra > 80]

print(Dor30_SFR)


plt.plot(Dor30_RP, Dor30_SFR, 'bo')
plt.title("RP vs SFR (rotation included)")
plt.xlabel("Ram Pressure",fontsize=16)
plt.ylabel("SFR",fontsize=16)
plt.tight_layout()
plt.savefig("SFR_RP_30Dor.pdf")
plt.close()

# dec range is -65 to -72
# ra range is 70 to 90

numbin = 4
raGrid = np.arange(70,90,(20/numbin))
decGrid = np.arange(-65,-72,-(7/numbin))

print("Testing SFR: ")
test = np.where((dec > -70) & (dec < -64) & (ra < 90) & (ra > 80))
print(SFR[test])


#rows, cols = (len(raGrid), len(decGrid))
#SFR_new = [[0 for i in range(cols)] for j in range(rows)]
SFR_new = [0 for k in range(len(raGrid)*len(decGrid))]
RP_new = [0 for k in range(len(raGrid)*len(decGrid))]
counter = 0
for i in range(0,len(raGrid)-1):
  for j in range(0,len(decGrid)-1):
      counter +=1
      with warnings.catch_warnings():
          warnings.simplefilter("ignore", category=RuntimeWarning)
          location = np.where((dec < decGrid[j]) & (dec > decGrid[j+1]) & (ra < raGrid[i+1]) & (ra > raGrid[i]))
          avgSFR = np.mean(10**(SFR[location]))
          avgRP = np.mean(RP[location])
          if (avgSFR > 0):
              SFR_new[counter] = np.log10(avgSFR)
              RP_new[counter] = np.log10(avgRP)
print("SFR averages: ")
print(SFR_new)

xnew = raGrid + (20/numbin)/2
ynew = decGrid -(7/numbin)/2

print(xnew)
print(ynew)

x,y = np.meshgrid(xnew, ynew)

plt.scatter(x,y, c = SFR_new,cmap='Blues_r',vmin = -8.0, vmax = -6)
plt.gca().invert_xaxis()
plt.xlabel("RA")
plt.ylabel("DEC")
plt.colorbar()
plt.title(r"Average SFR")
plt.savefig("AverageSFR_coarseRes.pdf")
plt.close()


plt.scatter(x,y, c = RP_new,cmap='Greens_r',vmin = -13, vmax = -12.5)
plt.gca().invert_xaxis()
plt.xlabel("RA")
plt.ylabel("DEC")
plt.colorbar()
plt.title(r"Average RP")
plt.savefig("AverageRP_coarseRes.pdf")
plt.close()

plt.plot(RP_new, SFR_new,'ko')
plt.xlim(-13,-12.5)
plt.ylim(-8,-6)
plt.xlabel("RP")
plt.ylabel("SFR")
plt.title(r"coarse-grained")
plt.savefig("RPSFR_coarseRes.pdf")
plt.close()


