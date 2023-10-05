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
#eddy = (3.0856e21*0.667/5e6)/3.155e13
eddy = 0.667/(0.5*0.11)

def _pcr(field, data):
  return 0.3333*data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")


yt.add_field(('gas', u'pcr'), function = _pcr, units="g/(cm*s**2)",display_name=r"P$_{cr}$", dimensions=dimensions.pressure)



plotvar = 'pcr'
#varmax = 1.e-24
varmax = 1.e-15
#varmax = 2.3e-27
varmin = 1.e-16

ts = yt.DatasetSeries('../cr.out1*',parallel=5)
for ds in ts.piter():
 print(ds.field_list)
 time = ds.current_time.v/eddy
# time = ds.current_time.v/eddy
 slc = yt.SlicePlot(ds, 'z', plotvar,fontsize=26)
 slc.set_zlim(plotvar, varmin, varmax)
 slc.set_cmap(field=plotvar, cmap='plasma')
 #slc.hide_colorbar()
 #slc.hide_axes()
 slc.set_xlabel('x')
 slc.set_ylabel('y')
# slc.set_width((2.0*kpc, 2.0*kpc))
# slc.set_log(plotvar, False)
# slc.annotate_streamlines('Bcc1','Bcc2')
# slc.annotate_quiver('vel1','vel2')
# slc.annotate_title("t = %3.0f Myr" % time)
 slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
# slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
 slc.save()
