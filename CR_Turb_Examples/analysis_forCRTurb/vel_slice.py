# plots slices of density for various times

import yt
from yt.units import dimensions
from yt.units import pc
from yt import YTQuantity
yt.enable_parallelism()

# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.
eddy = (3.0856e21/4e6)/3.155e13


def _velMag(field, data):
  return 1.E8*data['velocity_magnitude']

def _velx(field, data):
  return 1.E8*data['velocity_x']

yt.add_field(('gas', u'velx'), function = _velx, units="cm/s",display_name=r"v$_{x}$")
yt.add_field(('gas', u'velMag'), function = _velMag, units="cm/s",display_name=r"v$_{tot}$")



plotvar = 'velx'
#varmax = 1.e-24
#varmax = 5.e5
#varmax = 2.3e-27
#varmin = 5.e7

ts = yt.DatasetSeries('../cr.out1*',parallel=10)
for ds in ts.piter():
# time = ds.current_time.v/eddy
 time = ds.current_time
 ad = ds.all_data()
 print(ad.quantities.extrema("velx"))
 slc = yt.SlicePlot(ds, 'x', plotvar,fontsize=20)
# slc.set_zlim(plotvar, varmin, varmax)
 slc.set_cmap(field=plotvar, cmap='dusk')
 slc.set_xlabel('y/L')
 slc.set_ylabel('z/L')
# slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
 slc.annotate_title(r"t = %3.1f Myrs" % time)
 slc.save()
