import yt
import trident
from trident import LightRay
import aplpy
import numpy as np
#import matplotlib as plt
from yt.units.yt_array import YTQuantity
from yt import YTArray
import h5py
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import wcs
from astropy.wcs import WCS
import csv
import pandas as pd

inputArr = np.genfromtxt('ullyses_lmc_dr5_chadedits.csv',delimiter=',',invalid_raise=False)
ra_Yong = inputArr[:,0]
ra_Yong = ra_Yong[1:len(ra_Yong)]
dec_Yong = inputArr[:,1]
dec_Yong = dec_Yong[1:len(dec_Yong)]
#s_SFR = inputArr[:,2]
#s_SFR = s_SFR[1:len(s_SFR)]

f = fits.open("LMC_ram_RP_rotate_Jan2022_newSkyCenter.fits")
w = wcs.WCS(f[0].header)
print(w)

arr = f["ramPres"].data

ra_test, dec_test = w.wcs_world2pix(ra_Yong, dec_Yong, 1)
plt.plot(ra_Yong,ra_test)
plt.xlabel("RA Inputs")
plt.ylabel("Pixel Outputs")
plt.savefig("Ra_Testing.pdf")
plt.close()

plt.plot(dec_Yong,dec_test)
plt.xlabel("DEC Inputs")
plt.ylabel("Pixel Outputs")
plt.savefig("DEC_Testing.pdf")
plt.close()


#arr = np.flipud(arr)
pixel_x = []
pixel_y = []
#plt.imsave('test.png',arr,cmap="kamae_r")

fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(projection=w)
plt.imsave('test.png',arr, origin='lower', cmap='kamae_r')
plt.xlabel(r'RA')
plt.ylabel(r'Dec')

overlay = ax.get_coords_overlay('icrs')
overlay.grid(color='white', ls='dotted')

RP_YongCoords = []
for j in range(0,len(ra_Yong)):
  ra, dec = w.wcs_world2pix(ra_Yong[j], dec_Yong[j], 1)
 # ra, dec = w.world_to_pixel(ra_Yong[j], dec_Yong[j]) #, 1))
  ra = int(ra)
  dec = int(dec)
  print("Celestial Coordinate inputs: ")
  print(ra_Yong[j])
  print(dec_Yong[j])
  print("Pixel Coordinate outputs: ")
  print(ra)
  print(dec)
  print("RP Value: ")
 # print(arr[ra,512-dec])
 # RP_YongCoords.append(arr[ra,512-dec])
  print(arr[ra,dec])
  RP_YongCoords.append(arr[dec,ra])

RA, Dec = w.wcs_world2pix(ra_Yong, dec_Yong, 1)
newRA, newDec = w.wcs_pix2world(RA, Dec, 1)

plt.plot(ra_Yong,RA,'ko')
plt.savefig("conversionRA.pdf")
plt.close()

plt.plot(dec_Yong,Dec,'ko')
plt.savefig("conversionDEC.pdf")
plt.close()

df = pd.read_csv('ullyses_lmc_dr5_chadedits.csv')
df['RP(dyne/cm^2)'] = RP_YongCoords
print(df)

sc = plt.scatter(ra_Yong,dec_Yong,c = RP_YongCoords,cmap="Greens")
plt.colorbar(sc)
plt.title("Ram Pressure (Rotation Included)")
plt.gca().invert_xaxis()
plt.xlabel("RA (ICRS)",fontsize=16)
plt.ylabel("DEC (ICRS)",fontsize=16)
plt.tight_layout()
plt.savefig("Scatter_RP_Sept2022.pdf")
plt.close()


df.to_csv('RP_Rotation_dr5_Sept2022.csv')

"""
with open('LMC_SFR_surface_density_ullyses_04052021.csv','r') as csvinput:
    with open('RP.csv', 'w') as csvoutput:
        writer = csv.writer(csvoutput, lineterminator='\n')
        reader = csv.reader(csvinput)

        all = []
        row = next(reader)
        row.append('RP (dyne/cm^2)')
        all.append(row)

        for row in reader:
            row.append(RP_YongCoords[0])
            writer.writerows(row)
           # all.append(RP_YongCoords[0])
"""

