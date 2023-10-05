import matplotlib.pyplot as plt
import yt
from yt.units import erg, pc
from yt.units import dimensions
from yt.units.yt_array import YTQuantity
import numpy as np
# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.

yt.enable_parallelism()

base="cr.out1."
 
 
times = []
totalDampArr = []  # array holding the streaming loss
ts = yt.DatasetSeries('../cr.out1*',parallel=False)
i = 0
for ds in ts.piter():
  ad = ds.all_data()
  time = ds.current_time.v/Myr                       # store the time (with units removed)
  times.append(time)
 
  totalDamp = ad.quantities.total_quantity("alfven_speed")*1.E8
  averageAlfven = ad.quantities.weighted_average_quantity("alfven_speed",weight = "ones")*1.E8
  totalDampArr.append(totalDamp)
  print(averageAlfven)
 
print("total Alfven speed")
print(totalDampArr)
 
plt.plot(times,totalDampArr[0:len(times)],label="Total Alfven Speed")
plt.xlabel("Time (Myrs)")
plt.ylabel(r"Alfven speed (cm s$^{-1}$)")
#plt.ylim(1.e-13,3.e-11)
plt.legend()
plt.tight_layout()
plt.savefig('AlfvenSpeed.png')
plt.close()

