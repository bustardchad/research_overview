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

eddy = 0.667/(0.5*0.11)

yt.enable_parallelism()

def _pcr(field, data):
 return ds.arr(data['Ec'].v*edenstocgs*0.3333, 'g/(cm*s**2)')

yt.add_field(('gas', u'pcr'), function = _pcr, units="g/(cm*s**2)", dimensions=dimensions.pressure)


 
times = []
totalDampArr = []  # array holding the streaming loss
ts = yt.DatasetSeries('../cr.out1.*',parallel=5)
i = 0
for ds in ts.piter():
  grad_fields = ds.add_gradient_fields(("gas","pcr"))

  # We don't need to do the same for the pressure field because yt already
  # has pressure gradient fields. Now, define the "degree of hydrostatic
  # equilibrium" field.


  def _misalign(field, data):
    # Remember that g is the negative of the potential gradient
    gx = data[("gas", "pcr_gradient_x")]
    bx = data[("gas", "magnetic_field_x")]
    gy = data[("gas", "pcr_gradient_y")]
    by = data[("gas", "magnetic_field_y")]
    gz = data[("gas", "pcr_gradient_z")]
    bz = data[("gas", "magnetic_field_z")]
    noalign = np.sqrt(gx * gx + gy * gy + gz * gz)*data['magnetic_field_magnitude']
    align = np.sqrt((gx*bx)**2.0 + (gy*by)**2.0 + (gz*bz)**2.0)
    h = align/noalign
    return h


  ds.add_field(
    ("gas", "misalign"),
    function=_misalign,
    units="",
    take_log=True,
    display_name=r"$\frac{|v_{A} \cdot \nabla P_{CR}|}{|v_{A} \nabla P_{CR}|}$",
    sampling_type="cell",
  )



  # The gradient operator requires periodic boundaries.  This dataset has
  # open boundary conditions.
 # ds.force_periodicity()
  ad = ds.all_data()
  time = ds.current_time.v/eddy                       # store the time (with units removed)
  times.append(time)
  
  slc = yt.ProjectionPlot(ds, "z", "misalign",weight_field="ones",fontsize=26)
  slc.set_zlim("misalign",1e-2,1e0)
  slc.set_cmap(field="misalign", cmap='kamae')
  slc.set_xlabel('x')
  slc.set_ylabel('y')
  slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
  slc.save()
 
