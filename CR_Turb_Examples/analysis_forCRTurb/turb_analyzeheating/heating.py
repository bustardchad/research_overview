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
 
yt.enable_parallelism()
 
# Collisional energy losses
def collisions(field,data):
          return np.abs(data['user_out_var2'])*YTQuantity(1,"erg/cm**3/s")
 
 
yt.add_field(("gas","collisions"), function=collisions,display_name="Collisional Energy Loss Rate",units="erg/cm**3/s")

#Streaming energy losses
def heating(field,data):
          return np.abs(data['user_out_var7'])*(edenstocgs/3.155e13)*YTQuantity(1,"erg/cm**3/s")
 
 
yt.add_field(("gas","heating"), function=heating,display_name="Streaming Energy Loss Rate",units="erg/cm**3/s") 
 
times = []
totalCollArr = []  # array holding the collisional energy loss
totalDampArr = []  # array holding the streaming loss
totalCRLossArr = [] # array holding the total (collisional + streaming) loss
ts = yt.DatasetSeries('../cr.out2*',parallel=False)
i = 0

for ds in ts.piter():
  ad = ds.all_data()
  time = ds.current_time.v/Myr                       # store the time (with units removed)
  times.append(time)
 
  totalDamp = ad.quantities.total_quantity("heating")
 
  totalDampArr.append(totalDamp)
 
print("Collisionless")
print(totalDampArr)
 
plt.plot(times,totalDampArr[0:len(times)],label="Collisionless")
plt.xlabel("Time (Myrs)")
plt.ylabel(r"CR Energy Loss Rate (erg cm$^{-3}$ s$^{-1}$)")
#plt.ylim(1.e-13,3.e-11)
plt.legend()
plt.tight_layout()
plt.savefig('TotalCRLoss.png')
plt.close()

