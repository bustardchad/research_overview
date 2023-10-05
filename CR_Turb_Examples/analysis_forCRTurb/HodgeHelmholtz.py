# Hodge-Helmholtz decomposition of 3D velocity field into compressive and solenoidal components
# 
# Outline:
#	1. v(r) --> v(k) by taking a 3D FFT of velocity field v(r)
#	2. Find solenoidal component (k \cdot vsol(k) = 0) as:
#		For each i in 1 to 3:
#			visol(k) = sum(delta_ij - kikj/k^2)vj(k) for j = 1 to 3
#
# 	3. Find compressive component as vicomp(k) = vi(k) - visol(k)
#	4. Produce a power spectrum of each component, compare integrated power, etc.
#	5. (Optionally) Project back into physical space via inverse FFT


import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units.yt_array import YTQuantity
from yt.fields.api import ValidateParameter


# import helper functions
from utils.spectra import fft_comp as fft_comp
from utils.spectra import fft_comp_power as fft_comp_power
from utils.HH_utils import decompose, create_spectra


# TODO make shape of datacube, file_path, and files arguments passed in at command line
file_path = "files/"
files = "cr.out1.0004*"

ts = yt.DatasetSeries(file_path+files)
vtot_all_times = np.zeros(255)
vsol_all_times = np.zeros(255)
vcomp_all_times = np.zeros(255)

for ds in ts:
    dd = ds.all_data()

    # Decompose and return k and spectra for v_tot, v_sol, v_comp
    vx, vy, vz, v_comp_x, v_comp_y, v_comp_z, v_sol1_x, v_sol1_y, v_sol1_z = decompose(ds)

    v_sol_x = vx - v_comp_x
    v_sol_y = vy - v_comp_y
    v_sol_z = vz - v_comp_z


    power_ratio = (v_sol_x.var() + v_sol_y.var() + v_sol_z.var())/(v_comp_x.var() + v_comp_y.var() + v_comp_z.var())
    #power_ratio = (v_sol_x.mean()**2.0 + v_sol_y.mean() + v_sol_z.mean()**2.0)/(v_comp_x.mean()**2.0 + v_comp_y.mean()**2.0 + v_comp_z.mean()**2.0)
    
    print("Variance in each component")
    print('Solenoidal x, y, z: ' + str(v_sol_x.var()) + ", " + str(v_sol_y.var()) + ", " + str(v_sol_z.var()))
    print('Compressive x, y, z: ' + str(v_comp_x.var()) + ", " + str(v_comp_y.var()) + ", " + str(v_comp_z.var()))
    
    print("Average squared in each component")
    print('Solenoidal x, y, z: ' + str(v_sol_x.mean()**2.0) + ", " + str(v_sol_y.mean()**2.0) + ", " + str(v_sol_z.mean()**2.0))
    print('Compressive x, y, z: ' + str(v_comp_x.mean()**2.0) + ", " + str(v_comp_y.mean()**2.0) + ", " + str(v_comp_z.mean()**2.0))
 
    print("Ratio of solenoidal to compressive power: " + str(power_ratio))  

    # Calculate power spectra of total, solenoidal, and compressive kinetic energy components
    kvec, tot_KE_spec, sol_KE_spec, comp_KE_spec = create_spectra(ds, v_sol_x, v_sol_y, v_sol_z, v_comp_x, v_comp_y, v_comp_z)

    # Stack together different time snapshots to later do time-series analysis
    vtot_all_times = np.vstack([vtot_all_times, tot_KE_spec])
    vsol_all_times = np.vstack([vsol_all_times, sol_KE_spec])
    vcomp_all_times = np.vstack([vcomp_all_times, comp_KE_spec])


# Space for time-series analysis -- i.e. averaging over many snapshots, etc.

#vtot_all_times = vtot_all_times[1:8,:]  # last few rows (times) of spectrum array
#vsol_all_times = vsol_all_times[1:8,:]  # last few rows (times) of spectrum array
#vcomp_all_times = vcomp_all_times[1:8,:]  # last few rows (times) of spectrum array

# take mean, min, max over time (axis = 0)
avg_tot = np.mean(vtot_all_times,axis = 0)
mina_tot = np.min(vtot_all_times,axis = 0)
maxa_tot = np.max(vtot_all_times,axis = 0)

avg_sol = np.mean(vsol_all_times,axis = 0)
mina_sol = np.min(vsol_all_times,axis = 0)
maxa_sol = np.max(vsol_all_times,axis = 0)

avg_comp = np.mean(vcomp_all_times,axis = 0)
mina_comp = np.min(vcomp_all_times,axis = 0)
maxa_comp = np.max(vcomp_all_times,axis = 0)

k = kvec[0:len(avg_tot)]

print("Plotting average spectra (with min and max) over several snapshots: ")

# normalize spectra by average value of total velocity spectrum at 3rd mode
norma = avg_tot[2]*k[2]**2

# Plot spectra together
# (Optional) Uncomment "fill_between" lines to show a shaded region between min and max power for each line
plt.loglog((k), avg_tot*(k**2)/norma, 'ko-',label=r"Total (Compressive + Solenoidal)")
#plt.fill_between((k), mina_tot*(k**2), maxa_tot*(k**2), facecolor='black', alpha=0.5)
plt.loglog((k), avg_sol*(k**2)/norma, 'bo-',label=r"Solenoidal")
#plt.fill_between((k), mina_sol*(k**2), maxa_sol*(k**2), facecolor='blue', alpha=0.5)
plt.loglog((k), avg_comp*(k**2)/norma, 'go-',label=r"Compressive")
#plt.fill_between((k), mina_comp*(k**2), maxa_comp*(k**2), facecolor='green', alpha=0.5)
plt.ylim(1E-4,5E0)
#plt.xlim(3E-2,2E0)
plt.xlim(0.5,50)
plt.xlabel(r"$k$",fontsize=18)
plt.ylabel(r"Power Spectra x k^2",fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig("HH_figure.pdf")
plt.close()

print("Total Values: ")
print(avg_tot)
print(mina_tot)
print(maxa_tot)
print("Solenoidal Values: ")
print(avg_sol)
print(mina_sol)
print(maxa_sol)
print("Compressive Values: ")
print(avg_comp)
print(mina_comp)
print(maxa_comp)
print("k: ")
print(k)

