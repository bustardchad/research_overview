# Calculates broad turbulent fluctuation statistics like delta rho / rho, delta Ec/ Ec, etc.

import yt
from yt.units import dimensions
from yt.units import pc
from yt import YTQuantity
import numpy as np
yt.enable_parallelism()


devdensArr = []
devvelArr = []
devEcArr = []

file_path = "../"
files = "cr.out1.0004*"

ts = yt.DatasetSeries(file_path+files,parallel=False)

for ds in ts.piter():
  dd = ds.all_data()

  # Calculate average field quantities
  avg_density = dd.quantities.weighted_average_quantity("density",weight="ones")
  avg_vel = dd.quantities.weighted_average_quantity("velocity_magnitude",weight="ones")
  avg_Ec = dd.quantities.weighted_average_quantity("Ec",weight="ones")

  # Since compressive turbulence leads to lognormal distributions (usually), 
  # we move this calculation into log space, i.e. we'll take the log10(density / average density)
  # and then calculate it's standard deviation

  def density_log(field, data):
    return np.log(data['density']/avg_density)

  ds.add_field(('gas', u'density_log'), function = density_log, units="")

  def vel_log(field, data):
    return np.log(data['velocity_magnitude']/avg_vel)

  ds.add_field(('gas', u'vel_log'), function = vel_log, units="")

  def Ec_log(field, data):
    return np.log(data['Ec']/avg_Ec)

  ds.add_field(('gas', u'Ec_log'), function = Ec_log, units="")

  # Take the standard deviation, which in yt lingo is very incorrectly called the "variance"
  stddev_density, avg_density = dd.quantities.weighted_variance("density_log",weight="ones")
  stddev_vel,avg_vel = dd.quantities.weighted_variance("vel_log",weight="ones")
  stddev_Ec,avg_Ec = dd.quantities.weighted_variance("Ec_log",weight="ones")
  devdensArr.append(stddev_density)
  devvelArr.append(stddev_vel)
  devEcArr.append(stddev_Ec)

  print("delta rho ")
  print(devdensArr)
  print("delta v ")
  print(devvelArr)
  print("delta Ec ")
  print(devEcArr)
