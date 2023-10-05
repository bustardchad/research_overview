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

def _densratio(field, data):
  return data['density']/(YTQuantity(1,"g/cm**3")*0.03066)


yt.add_field(('gas', u'densratio'), function = _densratio, units="",display_name=r"$\rho / \rho_{i}$", dimensions=dimensions.density)



plotvar = 'densratio'
#varmax = 1.e-24
varmax = 3.0
#varmax = 2.3e-27
varmin = 0.3

ts = yt.DatasetSeries('../cr.out1*',parallel=5)
for ds in ts.piter():
 dd = ds.all_data()
 print(ds.field_list)
 time = ds.current_time.v/eddy
# time = ds.current_time.v/eddy
 print(dd.quantities.extrema("densratio"))
 slc = yt.SlicePlot(ds, 'z', plotvar,fontsize=26)
 slc.set_zlim(plotvar, varmin, varmax)
 slc.set_log(plotvar,False)
 slc.set_cmap(field=plotvar, cmap='kelp')
 #slc.hide_colorbar()
 slc.hide_axes()
 slc.set_xlabel('x/L')
 slc.set_ylabel('y/L')
# slc.set_width((2.0*kpc, 2.0*kpc))
# slc.set_log(plotvar, False)
# slc.annotate_streamlines('Bcc1','Bcc2')
# slc.annotate_quiver('vel1','vel2')
# slc.annotate_title("t = %3.0f Myr" % time)
 slc.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), kernellen = 100, lim=(0.5, 0.65),alpha = 0.6)
 slc.annotate_title(r"t = %3.1f $\tau_{\rm eddy}$" % time)
# plot = slc.plots[('gas',plotvar)]
# colorbar = plot.cb
# slc._setup_plots()
# colorbar.set_ticks([0.3, 1.0, 3.0])
# colorbar.set_ticklabels(['0.3', '1.0', '3.0'])
 #slc.set_cbar_minorticks(plotvar, "on")
# slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
 slc.save(suffix='pdf')
