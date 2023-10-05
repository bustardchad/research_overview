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


 
times = []
totalDampArr = []  # array holding the streaming loss
ts = yt.DatasetSeries('../cr.out1.*',parallel=False)
i = 0
for ds in ts.piter():
  grad_fields = ds.add_gradient_fields(("gas","pcr"))

  # We don't need to do the same for the pressure field because yt already
  # has pressure gradient fields. Now, define the "degree of hydrostatic
  # equilibrium" field.


  def _absGradPcr(field, data):
    # Remember that g is the negative of the potential gradient
    gx = data[("gas", "pcr_gradient_x")]
    gy = data[("gas", "pcr_gradient_y")]
    gz = data[("gas", "pcr_gradient_z")]
    h = np.sqrt(gx * gx + gy * gy + gz * gz)*(1.0/3.0856e21)
    return h


  ds.add_field(
    ("gas", "absGradPcr"),
    function=_absGradPcr,
    units="g/(cm**2 * s**2)",
    take_log=False,
    display_name="Grad Pcr",
    sampling_type="cell",
  )

  # The gradient operator requires periodic boundaries.  This dataset has
  # open boundary conditions.
 # ds.force_periodicity()
  ad = ds.all_data()
  time = ds.current_time.v/Myr                       # store the time (with units removed)
  times.append(time)
 
  totalDamp = ad.quantities.total_quantity("absGradPcr")
 
  totalDampArr.append(totalDamp)
 
 
print("total Grad Pcr")
print(totalDampArr)
 
plt.plot(times,totalDampArr[0:len(times)],label="Total Grad Pcr")
plt.xlabel("Time (Myrs)")
plt.ylabel(r"|$\nabla P_{CR}$|")
#plt.ylim(1.e-13,3.e-11)
plt.legend()
plt.tight_layout()
plt.savefig('TotalGradPcr.png')
plt.close()

