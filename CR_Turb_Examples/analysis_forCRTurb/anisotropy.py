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

from utils.spectra import fft_comp as fft_comp
from utils.spectra import fft_comp_power as fft_comp_power

def velocity_spectrum(ds):
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
    v_k_x = fft_comp_power(ds, ("gas","velocity_x"), max_level, low, dims)
    v_k_y = fft_comp_power(ds, ("gas","velocity_y"), max_level, low, dims)
    v_k_z = fft_comp_power(ds, ("gas","velocity_z"), max_level, low, dims)
 

    print("Shapes of v_k componennts: ")
    print(v_k_x.shape)   
    print(v_k_y.shape)   
    print(v_k_z.shape)
   
    # wavenumbers
    L = (ds.domain_right_edge - ds.domain_left_edge).d
    #kx = np.fft.rfftfreq(nx).reshape(nx,1,1)
    #ky = np.fft.rfftfreq(ny).reshape(ny,1,1)
    #kz = np.fft.rfftfreq(nz)
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
    #k2 = kx3d**2 + ky3d**2 + kz3d**2
    k2 = kx3d**2 + ky3d**2 + kz3d**2
    k = np.sqrt(k2)
 
    whichbin = np.digitize(k.flat, kbins)
    ncount = np.bincount(whichbin)
    
    v_x_spectrum = np.zeros(len(ncount)-1)
    v_y_spectrum = np.zeros(len(ncount)-1)
    v_z_spectrum = np.zeros(len(ncount)-1)

    for n in range(0,len(ncount)-1):
        v_x_spectrum[n] = np.sum(v_k_x.flat[whichbin==n])
        v_y_spectrum[n] = np.sum(v_k_y.flat[whichbin==n])
        v_z_spectrum[n] = np.sum(v_k_z.flat[whichbin==n])
    
    k = kbins[0:N]
    v_x_spectrum = v_x_spectrum[0:N]
    v_y_spectrum = v_y_spectrum[0:N]
    v_z_spectrum = v_z_spectrum[0:N]
    
    return k, v_x_spectrum, v_y_spectrum, v_z_spectrum



def doit(ds,v_sol_x, v_sol_y, v_sol_z, v_comp_x, v_comp_y, v_comp_z):

    # a FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this.

    max_level = ds.index.max_level
    print(max_level)
 #   ref = int(np.product(ds.ref_factors[0:max_level]))

    low = ds.domain_left_edge
    print(low)
   # dims = ds.domain_dimensions*ref
    dims = ds.domain_dimensions*int(1.0)
    print(dims)
    nx, ny, nz = dims

    nindex_rho = 1./2.
   # nindex_rho = 0.0

    #Kk = np.zeros( (nx//2+1,ny//2+1) )
    Kk = np.zeros( (nx//2+1, ny//2+1, nz//2+1))
    Kk_sol = np.zeros( (nx//2+1, ny//2+1, nz//2+1))
    Kk_comp = np.zeros( (nx//2+1, ny//2+1, nz//2+1))
    
    for vel in [("gas","velocity_x"), ("gas","velocity_y"), ("gas","velocity_z")]:
    	Kk += fft_comp_power(ds, vel, max_level, low, dims)
    
    for vel in [v_sol_x, v_sol_y, v_sol_z]:
    	ru = np.fft.fftn(vel)[0:nx//2+1,0:ny//2+1,0:nz//2+1]
    	ru = 8.0*ru/(nx*ny*nz)
    	Kk_sol += np.abs(ru)**2.0
    
    for vel in [v_comp_x, v_comp_y, v_comp_z]:
    	ru = np.fft.fftn(vel)[0:nx//2+1,0:ny//2+1,0:nz//2+1]
    	# ru = np.fft.rfft(rho**nindex_rho * u)[0:nx//2+1,0:ny//2+1,0]
    	ru = 8.0*ru/(nx*ny*nz)
    	Kk_comp += np.abs(ru)**2.0


    # wavenumbers
    L = (ds.domain_right_edge - ds.domain_left_edge).d
    print(L)
    print(np.fft.rfftfreq(nx))
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



ts = yt.DatasetSeries("../cr.out1.0004*")
vx_all_times = np.zeros(255)
vy_all_times = np.zeros(255)
vz_all_times = np.zeros(255)


for ds in ts:
    dd = ds.all_data()
    
    # Decompose and return k and spectra for v_tot, v_sol, v_comp
    kvec, Kk_vx, Kk_vy, Kk_vz = velocity_spectrum(ds)

    vx_all_times = np.vstack([vx_all_times, Kk_vx])
    vy_all_times = np.vstack([vy_all_times, Kk_vy])
    vz_all_times = np.vstack([vz_all_times, Kk_vz])

# take mean, min, max over time (axis = 0)
avg_x = np.mean(vx_all_times,axis = 0)
mina_x = np.min(vx_all_times,axis = 0)
maxa_x = np.max(vx_all_times,axis = 0)

avg_y = np.mean(vy_all_times,axis = 0)
mina_y = np.min(vy_all_times,axis = 0)
maxa_y = np.max(vy_all_times,axis = 0)

avg_z = np.mean(vz_all_times,axis = 0)
mina_z = np.min(vz_all_times,axis = 0)
maxa_z = np.max(vz_all_times,axis = 0)

k = kvec[0:len(avg_x)]


print("V_x Values: ")
print(avg_x)
print(mina_x)
print(maxa_x)
print("V_y Values: ")
print(avg_y)
print(mina_y)
print(maxa_y)
print("V_z Values: ")
print(avg_z)
print(mina_z)
print(maxa_z)
print("V_x/V_tot Values: ")
print(avg_x/(avg_x+avg_y+avg_z))
print("k: ")
print(k)
