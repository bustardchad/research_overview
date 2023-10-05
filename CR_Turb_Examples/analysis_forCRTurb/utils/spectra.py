import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units.yt_array import YTQuantity

# Calculates the Fourier transform of an input field "iu"
def fft_comp(ds, iu, level, low, delta ):

    # FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this in general, which allows this script to run on 
    # simulation grids with non-uniform resolution. In this specific application,
    # our grid is already uniform so this next line is overkill.
    cube = ds.covering_grid(level, left_edge=low,
                            dims=delta,
                            fields=[iu])

    u = cube[iu].d

    ru = np.fft.fftn(u) # computes 3D FFT over all axes
    return u, ru # iu, fourier component of iu 

def fft_comp_power(ds, iu, level, low, delta ):

    # FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this in general, which allows this script to run on 
    # simulation grids with non-uniform resolution. In this specific application,
    # our grid is already uniform so this next line is overkill.
    cube = ds.covering_grid(level, left_edge=low,
                            dims=delta,
                            fields=[iu])

    u = cube[iu].d

    nx, ny, nz = u.shape

    # We keep just the first half of the axes in each direction.  Our
    # normalization has an '8' to account for this clipping to one
    # octant.
    ru = np.fft.fftn(u)[0:nx//2+1,0:ny//2+1,0:nz//2+1]
    ru = 8.0*ru/(nx*ny*nz)

    return np.abs(ru)**2  # absolute value and square to get the power

# Full version that takes in multiple fields (irho, iu) and additional input nindex_rho
def fft_comp_withrho(ds, irho, iu, nindex_rho, level, low, delta ):

    # FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this in general, which allows this script to run on 
    # simulation grids with non-uniform resolution. In this specific application,
    # our grid is already uniform so this next line is overkill.
    cube = ds.covering_grid(level, left_edge=low,
                            dims=delta,
                            fields=[irho, iu])

    rho = cube[irho].d
    u = cube[iu].d

    nx, ny, nz = rho.shape

    # Take FFT of rho^nindex_rho x u -- for kinetic energy, nindex_rho should = 1/2
    ru = np.fft.fftn(rho**nindex_rho * u)[0:nx//2+1,0:ny//2+1,0:nz//2+1]
    ru = 8.0*ru/(nx*ny*nz)

    return np.abs(ru)**2  # e.g. gives rho v^2 for kinetic energy