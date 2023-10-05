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

def _pcr(field, data):
 return ds.arr(data['Ec'].v*edenstocgs*0.3333, 'g/(cm*s**2)')
 
yt.add_field(('gas', u'pcr'), function = _pcr, units="g/(cm*s**2)", dimensions=dimensions.pressure)
 
base="cr.out1."
 
plotvar = 'pcr'

 
 
times = []
totalDampArr = []  # array holding the streaming loss
ts = yt.DatasetSeries('../cr.out1*',parallel=False)
i = 0
for ds in ts.piter():
  ad = ds.all_data()
  time = ds.current_time.v/Myr                       # store the time (with units removed)
  times.append(time)
 
  totalDamp = ad.quantities.total_quantity("pcr")
 
  totalDampArr.append(totalDamp)
 
 
print("total Pcr")
print(totalDampArr)
 
plt.plot(times,totalDampArr[0:len(times)],label="Total Pcr")
plt.xlabel("Time (Myrs)")
plt.ylabel(r"CR Pressure (g cm$^{-1}$ s$^{-2}$)")
#plt.ylim(1.e-13,3.e-11)
plt.legend()
plt.tight_layout()
plt.savefig('TotalCRPressure.png')
plt.close()

