import yt
from yt.units.yt_array import YTQuantity
from yt.fields.api import ValidateParameter
import numpy as np
import matplotlib as plt
from yt.units import Mpc
import collections
import glob

yt.enable_parallelism()

pUnit = YTQuantity(1, 'cm**2/s**2')

def CoolingmR(field,data):
  return -data['cloo']*YTQuantity(1,'erg/cm**3/s')/(4.0*np.pi*2.41e-10) # in mR

yt.add_field(("gas","CoolingmR"), function=CoolingmR,units="erg/cm**3/s")


# By using wildcards such as ? and * with the load command, we can load up a
# Time Series containing all of these datasets simultaneously.




num_procs = 1
fns = glob.glob('../Stream_hdf5_plt_cnt_**9')
fns.sort()
print(fns.sort())
#ts = yt.load('../Stream_hdf5_plt_cnt_????')

# Calculate and store density extrema for all datasets along with redshift
# in a data dictionary with entries as tuples

# Create an empty dictionary
data = {}

# Iterate through each dataset in the Time Series (using piter allows it
# to happen in parallel automatically across available processors)
for sto,fn in yt.parallel_objects(fns,num_procs,storage=data):
    ds = yt.load(fn)
    ad = ds.all_data()
    ds.periodicity = (True,True,True)
    ds.coordinates.x_axis[1] = 0
    ds.coordinates.y_axis[1] = 2
    dense_ad = ad.cut_region(['obj["CoolingmR"] > 0.0'])
    # c = yt.ProjectionPlot(ds, 1, "cloo",weight_field="density",fontsize=28)
    proj = ds.proj("CoolingmR",1,data_source=dense_ad)
    frb = proj.to_frb(width=(100,'kpc'),resolution=[256,256])
    my_image = np.array(frb["CoolingmR"])
 #   print(my_image)
    #filtered_image = np.array(filter(lambda num: num > 0.0, my_image))
    filtered_image = my_image[my_image > 0.0]
 #   print(filtered_image)
    maxCol = np.amax(filtered_image)
    avgCol = np.average(filtered_image)
    stdCol = np.std(filtered_image)
    #proj.save_as_dataset()

    # reload data set and print out some stuff
    #extrema = dense_ad.quantities.extrema('CoolingmR')
    total = dense_ad.quantities.total_quantity('CoolingmR')
    var = dense_ad.quantities.weighted_variance('CoolingmR',weight='ones')

    # Fill the dictionary with extrema and redshift information for each dataset
    time = ds.current_time.in_units('Myr')
    sto.result_id = ds.basename
    sto.result = (time, maxCol, avgCol, stdCol, total)
   # data[ds.basename] = (time, extrema, var)
# Convert dictionary to ordered dictionary to get the right order

"""
od = collections.OrderedDict(sorted(data.items()))

# Print out all the values we calculated.
print("Dataset      Time       Halpha Min      Halpha Max      Halpha variance        Halpha averaage")
print("---------------------------------------------------------")
for key, val in od.items():
    print("%s       %05.3f      %5.3g mR         %5.3g mR         %5.3g mR             %5.3g mR" % \
           (key, val[0], val[1][0], val[1][1], val[2][0], val[2][1]))
"""

file1 = open("CoolFile_256.txt","w")#write mode
#if yt.is_root()
file1.write("Dataset                          Time       Lum Max    Lum average    Lum std dev      Lum Total      \n")
file1.write("---------------------------------------------------------------------------------------------------------------------------------------    \n ")
for fn, val in sorted(data.items()):
    file1.write("%s       %05.3f      %5.3g mR     %5.3g mR     %5.3g mR         %5.3e      \n" % \
           (fn, val[0], val[1], val[2], val[3], val[4]))    

file1.close()
