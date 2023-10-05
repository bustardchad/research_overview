# plots slices of density for various times

import yt
from yt.units import dimensions
from yt.units import pc
from yt import YTQuantity
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
yt.enable_parallelism()

# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.
eddy = (3.0856e21/4e6)/3.155e13

vm1 = 40.0
vm80 = 80.0

ds1 = yt.load('../cr.out1.00001.athdf')
ds80 = yt.load('../../vm80/cr.out1.00001.athdf')

# Example for a 512 by 512 slice in the xy plane.
ds1r = ds1.r[::64j, ::64j, 0]
ds80r = ds80.r[::64j, ::64j, 0]

flux1 = np.sqrt(np.array(ds1r['Fc1'])**2.0 + np.array(ds1r['Fc2'])**2.0 + np.array(ds1r['Fc3'])**2.0)*vm1

flux80 = np.sqrt(np.array(ds80r['Fc1'])**2.0 + np.array(ds80r['Fc2'])**2.0 + np.array(ds80r['Fc3'])**2.0)*vm80


# Plot the difference
fig, ax = plt.subplots()
im = ax.imshow(np.array(flux1)/np.array(flux80), origin = 'lower',cmap='RdBu_r',vmin=0.1, vmax=10.0,norm=colors.LogNorm())
cbar = fig.colorbar(im, ax=ax, extend='both')
plt.savefig('fluxratio_80_40.png')
plt.close()



dens1 = np.array(ds1r['density'])
dens80 = np.array(ds80r['density'])

# Plot the difference
fig, ax = plt.subplots()
im = ax.imshow(np.array(dens1), origin = 'lower',cmap='RdBu_r',norm=colors.LogNorm())
cbar = fig.colorbar(im, ax=ax, extend='both')
plt.savefig('dens_vm1.png')
plt.close()

