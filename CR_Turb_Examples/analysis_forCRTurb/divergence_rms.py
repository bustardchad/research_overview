# plots slices of density for various times

import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units import dimensions
from yt.units import pc
from yt import YTQuantity
yt.enable_parallelism()

# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.
eddy = (3.0856e21/5e6)/3.155e13


def _velMag(field, data):
  return ((1.E8/3.0856e21)*(data['velocity_divergence']))**2.0 # normalized to u/L

yt.add_field(('gas', u'velMag'), function = _velMag, units="1/s**2",display_name=r"|$\nabla \cdot v$|")



plotvar = 'velMag'
#varmax = 1.e-24
varmax = 5.e7
#varmax = 2.3e-27
varmin = 5.e5
arr = []
times = []

ts = yt.DatasetSeries('../cr.out1*',parallel=10)
for ds in ts.piter():
  dd = ds.all_data()
  time = ds.current_time/Myr
  times.append(time)
  arr.append(np.sqrt(dd.quantities.weighted_average_quantity('velMag', 'ones'))*(3.0856e21*0.67/1e7)/(3.0))

print(arr)

plt.semilogy(times,arr)
plt.xlabel("Time (Myrs)")
plt.ylabel(r"($\nabla \cdot v$)$\tau_{w}$/3 ")
plt.tight_layout()
plt.savefig("velDiv_RMS_tw_over3.pdf")
plt.close()
