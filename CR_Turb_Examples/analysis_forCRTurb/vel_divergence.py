# Script that calculates and plots velocity divergence
# Also calculates average velocity divergence

import yt
from yt.units import dimensions
from yt.units import pc
from yt import YTQuantity
import numpy as np
yt.enable_parallelism()

# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.
eddy = (3.0856e21/5e6)/3.155e13


def _velDivRMS(field, data):
  return np.sqrt(data['velocity_divergence']*data['velocity_divergence'])

yt.add_field(('gas', u'velDivRMS'), function = _velDivRMS, units="1/s",display_name=r"<$\nabla \cdot v$>")

avgArr = []

plotvar = 'velDivRMS'
#varmax = 1.e-24
varmax = 5.e7
#varmax = 2.3e-27
varmin = 5.e5

ts = yt.DatasetSeries('../cr.out1.000*',parallel=False)
for ds in ts.piter():
# time = ds.current_time.v/eddy
 time = ds.current_time
 dd = ds.all_data()
 avgArr.append(dd.quantities.weighted_average_quantity("velDivRMS",weight="ones"))
 print(dd.quantities.weighted_average_quantity("velocity_magnitude",weight="ones"))
 slc = yt.SlicePlot(ds, 'z', plotvar,fontsize=20)
# slc.set_zlim(plotvar, varmin, varmax)
 slc.set_cmap(field=plotvar, cmap='dusk')
 slc.set_xlabel('x/L')
 slc.set_ylabel('y/L')
# slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
 slc.annotate_title(r"t = %3.1f Myrs" % time)
 slc.save()


print("Average RMS velocity divergence: ")
print(avgArr)
