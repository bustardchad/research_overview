import matplotlib.pyplot as plt
import yt
from yt.units import erg, pc
from yt.units.yt_array import YTQuantity
import numpy as np
# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.
eddy = 0.667/(0.5*0.11)
yt.enable_parallelism()

#Streaming energy losses
def crscale(field,data):
          return np.abs(data['user_out_var8'])/(0.667*3.0856e21)
 
 
yt.add_field(("gas","crscale"), function=crscale,display_name=r"$L_{CR}/L_{0}$",units="") 
 
times = []
totalCollArr = []  # array holding the collisional energy loss
totalDampArr = []  # array holding the streaming loss
totalDampDensArr = []  # array holding the streaming loss
totalCRLossArr = [] # array holding the total (collisional + streaming) loss
ts = yt.DatasetSeries('../cr.out2*',parallel=False)
i = 0

for ds in ts.piter():
  ad = ds.all_data()
  time = ds.current_time.v/eddy
  cut_ad = ad.cut_region(['obj["gas", "crscale"] < 1e29'])
  times.append(time)
 
  totalDamp = cut_ad.quantities.weighted_average_quantity("crscale",weight="ones")
 # totalDamp_dens = ad.quantities.weighted_average_quantity("crscale",weight="density")
  print("Average crscale: ")
  print(totalDamp)
  totalDampArr.append(totalDamp)
 # totalDampDensArr.append(totalDamp_dens)
   
 # slc = yt.ProjectionPlot(ds, "z", "crscale",weight_field="ones")
 # slc.set_zlim("crscale",1,1e3)
 # slc.set_cmap(field="crscale", cmap='hot_r')
 # slc.save()

  slc = yt.SlicePlot(ds, "z", "crscale",fontsize=26)
  slc.set_zlim("crscale",0.1,1e3)
  slc.set_cmap(field="crscale", cmap='hot_r')
  slc.set_xlabel('x')
  slc.set_ylabel('y')
  slc.hide_axes()
  slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
  slc.save(suffix='pdf')

print("lcr volume weighted average")
print(totalDampArr)

#print("lcr density weighted average")
#print(totalDampDensArr)

