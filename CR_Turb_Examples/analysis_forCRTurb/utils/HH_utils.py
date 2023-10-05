# Utilities for Hodge-Helmholtz decomposition
# 
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

from utils.spectra import fft_comp, fft_comp_power

def decompose(ds):
    # Calls FFT to do step 1
    # Does steps 2-4 for an individual snapshot and returns
    #	Vector of wavenumbers, total velocity power spectrum, solenoidal power spectrum, 
    #   compressive power spectrum
   

    # Step 1 ..........................................................................

    # a FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this.
    max_level = ds.index.max_level

    low = ds.domain_left_edge
    dims = ds.domain_dimensions*int(1.0)
    nx, ny, nz = dims

    
    # FFT of v_x, v_y, v_z
    vx, v_k_x = fft_comp(ds, ("gas","velocity_x"), max_level, low, dims)
    vy, v_k_y = fft_comp(ds, ("gas","velocity_y"), max_level, low, dims)
    vz, v_k_z = fft_comp(ds, ("gas","velocity_z"), max_level, low, dims)
 

    print("Shapes of v_k components: ")
    print(v_k_x.shape)   
    print(v_k_y.shape)   
    print(v_k_z.shape)
   
    # wavenumbers
    #kx = np.fft.rfftfreq(nx).reshape(nx,1,1)
    #ky = np.fft.rfftfreq(ny).reshape(ny,1,1)
    #kz = np.fft.rfftfreq(nz)
    kx = np.fft.fftfreq(nx)
    ky = np.fft.fftfreq(ny)
    kz = np.fft.fftfreq(nz)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
    k2 = kx3d**2 + ky3d**2 + kz3d**2
    k2[0,0,0] = 1. # get rid of infinite scale
   

    # new
    div_Vf_k = (v_k_x * kx3d + v_k_y * ky3d + v_k_z * kz3d)
    v_comp_overk = div_Vf_k / k2
    
    # Compressive components
    v_comp_x = np.fft.ifftn(v_comp_overk * kx3d)
    v_comp_y = np.fft.ifftn(v_comp_overk * ky3d)
    v_comp_z = np.fft.ifftn(v_comp_overk * kz3d)
    
    # Solenoidal components
    v_sol_x = np.fft.ifftn(v_k_x - (v_comp_overk * kx3d))
    v_sol_y = np.fft.ifftn(v_k_y - (v_comp_overk * ky3d))
    v_sol_z = np.fft.ifftn(v_k_z - (v_comp_overk * kz3d))

    return vx, vy, vz, v_comp_x, v_comp_y, v_comp_z, v_sol_x, v_sol_y, v_sol_z


def create_spectra(ds,v_sol_x, v_sol_y, v_sol_z, v_comp_x, v_comp_y, v_comp_z):

    # Only relevant if we have non-uniform resolution
    max_level = ds.index.max_level

    low = ds.domain_left_edge

    # Calculate dimensions of simulation domain in each direction
    dims = ds.domain_dimensions*int(1.0)
    nx, ny, nz = dims

    # Arrays to hold the Fourier power total, in solenoidal motions, 
    # and in compressive motions
    Kk = np.zeros( (nx//2+1, ny//2+1, nz//2+1))
    Kk_sol = np.zeros( (nx//2+1, ny//2+1, nz//2+1))
    Kk_comp = np.zeros( (nx//2+1, ny//2+1, nz//2+1))
    
    # Loop over velocities in x, y, z directions and calculate Fourier power
    for vel in [("gas","velocity_x"), ("gas","velocity_y"), ("gas","velocity_z")]:
    	Kk += fft_comp_power(ds, vel, max_level, low, dims)
    
    for vel in [v_sol_x, v_sol_y, v_sol_z]:
        ru = np.fft.fftn(vel)[0:nx//2+1,0:ny//2+1,0:nz//2+1]
        ru = 8.0*ru/(nx*ny*nz)
        Kk_sol += np.abs(ru)**2.0
    
    for vel in [v_comp_x, v_comp_y, v_comp_z]:
    	ru = np.fft.fftn(vel)[0:nx//2+1,0:ny//2+1,0:nz//2+1]
    	ru = 8.0*ru/(nx*ny*nz)
    	Kk_comp += np.abs(ru)**2.0


    # wavenumbers
    L = (ds.domain_right_edge - ds.domain_left_edge).d
    kx = np.fft.rfftfreq(nx)*nx/L[0]
    ky = np.fft.rfftfreq(ny)*ny/L[1]
    kz = np.fft.rfftfreq(nz)*nz/L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0/L)
    kmax = np.min(0.5*dims/L)

    kbins = np.arange(kmin, kmax, kmin)
    N = len(kbins)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

    whichbin = np.digitize(k.flat, kbins)
    ncount = np.bincount(whichbin)
    
    E_spectrum = np.zeros(len(ncount)-1)
    E_spectrum_sol = np.zeros(len(ncount)-1)
    E_spectrum_comp = np.zeros(len(ncount)-1)

    for n in range(0,len(ncount)-1):
        E_spectrum[n] = np.sum(Kk.flat[whichbin==n])
        E_spectrum_sol[n] = np.sum(Kk_sol.flat[whichbin==n])
        E_spectrum_comp[n] = np.sum(Kk_comp.flat[whichbin==n])

    k = kbins[0:N]
    E_spectrum = E_spectrum[0:N]
    E_spectrum_sol = E_spectrum_sol[0:N]
    E_spectrum_comp = E_spectrum_comp[0:N]


    return k, E_spectrum, E_spectrum_sol, E_spectrum_comp
