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
 
# compressions = v dot grad pc
def compression(field,data):
          return np.abs(data['user_out_var9'])*(edenstocgs/3.155e13)*YTQuantity(1,"erg/cm**3/s")


yt.add_field(("gas","compression"), function=compression,display_name=r"$v \cdot \nabla P_{CR}$",units="erg/cm**3/s")
 
yt.add_field(("gas","crscale"), function=crscale,display_name=r"$L_{CR}/L_{0}$",units="") 
 
times = []
totalCollArr = []  # array holding the collisional energy loss
totalDampArr = []  # array holding the streaming loss
totalDampDensArr = []  # array holding the streaming loss
totalCRLossArr = [] # array holding the total (collisional + streaming) loss
ts = yt.DatasetSeries('../cr.out2*',parallel=5)
i = 0

for ds in ts.piter():
  ad = ds.all_data()
  time = ds.current_time.v/eddy

  slc = yt.SlicePlot(ds, "z", "compression",fontsize=26)
  slc.set_zlim("compression",1e-30,1e-28)
  slc.set_cmap(field="compression", cmap='plasma')
  slc.set_xlabel('x')
  slc.set_ylabel('y')
  slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
  slc.save()


