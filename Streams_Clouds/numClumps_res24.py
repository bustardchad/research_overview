import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import cm

def add_legend(plot_args_lst, legend_loc = 'best', default_plot_args = {'color' : 'black'},
                  **kwargs):
       """Adds another legend to plot.
       Keyword Arguments:
       plot_args_lst      -- List of plotting arguments to show.
       legend_loc         -- Location of new legend (default 'best')
       default_plot_args  -- Arguments will be used for every item in new legend.
       kwargs             -- Will be passed to `plt.legend` of new legend.
       Example:
               add_legend([{'ls' : '-', 'label' : r'$\rho > \rho_{cl}/3$'}, {'ls' : '--', 'label' :   r'$T   < 2 T_{cl}$'}],
                                   default_plot_args = {'c' : 'k'},
                                   title = r'')

       Will add a legend with two different lines (both black).
       """
       ax = plt.gca()
       leg = ax.get_legend()
       linelst = []
       for cargs in plot_args_lst:
           for k, v in default_plot_args.items():
               if k not in cargs:
                   cargs[k] = v
           l, = plt.plot(np.nan, **cargs)
           linelst.append(l)
       o = kwargs.copy()
       if 'loc' not in o:
           o['loc'] = legend_loc
       if 'fontsize' not in o:
           o['fontsize'] = 7.5
       if legend_loc == 'above':
           o['loc'] = 'lower center'
           o['bbox_to_anchor'] = (0.5, 1.01)
       plt.legend(handles = linelst, **o)
       if leg is not None:
           ax.add_artist(leg) # Add old legend







Machs = [1,2,3,6]
paths_wind = ['Wind/res24/Mach1/plots/clumpDir/', 'Wind/res24/Mach2/plots/clumpDir/','Wind/res24/Mach3/plots/clumpDir/','Wind/res24/Mach6/plots/clumpDir/']

paths_shock = ['Shock/res24/Mach1/plots/clumpDir/', 'Shock/res24/Mach2/plots/clumpDir/','Shock/res24/Mach3/plots/clumpDir/', 'Shock/res24/Mach6/plots/clumpDir/']


# function that returns arrays for number of clumps and times
def clumps_vs_time(path):
    times = []
    numClumps = []
    # count the number of files in that path and loop over them in order
    count = 0
    for p in os.listdir(path):
        if os.path.isfile(os.path.join(path,p)):
            count +=1

    for i in range(0,10*count,10):
      print("File being read: " + str(path) + 'cloud.'+str(i).zfill(4)+'.hdf5')
      with pd.HDFStore(path+'cloud.'+str(i).zfill(4)+'.hdf5','r') as d: # read file in
        df = d.get("dat")
        dens = np.array(df[df.columns[11]])  # each clump has various properties, just need the sum of identified clumps...
        numClumps.append(len(dens)) # number of clumps
        times.append(0.1*i) # in units of cloud-crushing time
    return times, numClumps


colors = cm.plasma(np.linspace(0, 1, len(paths_wind)+1))
lw = 2

fig, axs = plt.subplots(2,1,sharex=True,sharey=False)

countWinds = 0
for path in paths_wind:
    times, numClumps = clumps_vs_time(path)
    axs[0].semilogy(times,numClumps,color=colors[countWinds],label=r"$\mathcal{M} = %1.1f $" %                Machs[countWinds],linewidth=lw,linestyle="dashed")
    countWinds +=1

countShocks = 0
for path in paths_shock:
    times, numClumps = clumps_vs_time(path)
    axs[1].semilogy(times,numClumps,color=colors[countShocks],linewidth=lw,linestyle="solid")
    countShocks +=1
    print(countShocks)

axs[1].set_xlabel(r"Time (t$_{cc}$)",fontsize=18)
axs[0].set_ylabel(r"# of Clumps",fontsize=18)
axs[1].set_ylabel(r"# of Clumps",fontsize=18)
axs[0].set_ylim(1,1e3)
axs[1].set_ylim(1,1e3)
axs[0].set_xlim(0,20)
axs[1].set_xlim(0,20)
axs[0].legend(loc='lower right')
add_legend([{'ls' : '-', 'label' : r'Shock'}, {'ls' : '--', 'label' :   r'Wind'}], default_plot_args = {'c' : 'k'},title = r'',loc='lower right')
plt.tight_layout()
plt.savefig("numClumps_chi100_Tcl1e4_varyMach_res24.pdf")
plt.close()

