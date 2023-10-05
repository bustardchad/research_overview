import yt
from yt.units.yt_array import YTQuantity
from yt.units import dimensions
import numpy as np
from numpy import *
import matplotlib as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import h5py
from yt.fields import interpolated_fields
from yt.fields.field_detector import \
    FieldDetector
from yt.utilities.linear_interpolators import \
    BilinearFieldInterpolator

yt.enable_parallelism()

pUnit = YTQuantity(1, 'cm**2/s**2')

# conversion factors
denstocgs = 6.85e-27
edenstocgs = 6.54e-11
Myr = 1.
kpc = 1.




def _d(field, data):
  return data['density']*denstocgs


yt.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)

profiles = []
labels = []
plot_specs = []

ds = yt.load('../cr.out1.00010.athdf')
dd = ds.all_data()
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profiles.append(yt.create_profile(dd, "d", ["cell_volume"],weight_field=None, fractional=True, accumulation=False,extrema={'cell_volume': (1e-3, 1e0),
                                      'd': (1e-29, 1e-27)}))

# Add labels
labels.append(r"$\beta \sim 1$")
plot_specs.append(dict(linestyle="-", linewidth=3))

ds = yt.load('../../higherBeta/diff3e26/cr.out1.00010.athdf')
dd = ds.all_data()
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profiles.append(yt.create_profile(dd, "d", ["cell_volume"],weight_field=None, fractional=True, accumulation=False,extrema={'cell_volume': (1e-3, 1e0),
                                      'd': (5e-29, 5e-27)}))

# Add labels
labels.append(r"$\beta \sim 100$")
plot_specs.append(dict(linestyle="--", linewidth=3))


################33
# Create the profile plot from the list of profiles.
plot = yt.ProfilePlot.from_profiles(profiles, labels=labels, plot_specs=plot_specs)
#plot.set_xlim(1.E-27,1.E-21)
#plot.set_ylim("ecr",1.E-3,1.E0)
#plot.set_ylim("gammaRay_lum",1.E-3,1.E0)
plot.set_ylabel("cell_volume", "Volume-Weighted Density PDF")
plot.annotate_title(r"Density Distributions")
#plot.annotate_title(r"L = 5, $\alpha$ = 1.5")
plot.save(suffix='pdf')
        
