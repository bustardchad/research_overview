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

def _ecr(field, data):
  return data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")

def _pcr_pg(field, data):  # pcr / rho cs^2
  return 0.333*data['Ec']*edenstocgs/((1e7**2.0)*data['density']*denstocgs) * YTQuantity(1,"g/cm**3")

yt.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)",display_name=r"E$_{cr}$", dimensions=dimensions.pressure)
yt.add_field(('gas', u'pcr_pg'), function = _pcr_pg, units="",display_name=r"P$_{cr}$/P$_{g}$")



plotvar = 'pcr_pg'
#varmax = 1.e-24
varmax = 5.e-2
#varmax = 2.3e-27
varmin = 5.e-3

ts = yt.DatasetSeries('../cr.out1*',parallel=False)
for ds in ts.piter():
 print(ds.field_list)
 time = ds.current_time.v/eddy
# time = ds.current_time.v/eddy
 slc = yt.ProjectionPlot(ds, 'x', plotvar,weight_field='ones',fontsize=26)
 slc.set_zlim(plotvar, varmin, varmax)
 slc.set_cmap(field=plotvar, cmap='dusk')
 #slc.hide_colorbar()
 #slc.hide_axes()
 slc.set_xlabel('y')
 slc.set_ylabel('z')
# slc.set_width((2.0*kpc, 2.0*kpc))
# slc.set_log(plotvar, False)
# slc.annotate_streamlines('Bcc1','Bcc2')
# slc.annotate_quiver('vel1','vel2')
# slc.annotate_title("t = %3.0f Myr" % time)
 slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
# slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
 slc.save()
