import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# at the same number of cloud-crushing times
num15 = []
time15 = [] # each step is 0.2t_cc
num3 = []
time3 = []  # each step is 0.4t_cc
num6 = []
time6 = []  # each step is 0.8t_cc

for i in range(0,60,10):
  Mach15 = 'res32/clumpDir/cloud.'+str(i).zfill(4)+'.hdf5'
  Mach3 = 'res64/clumpDir/cloud.'+str(i).zfill(4)+'.hdf5'
  Mach6 = 'res128/clumpDir/cloud.'+str(i).zfill(4)+'.hdf5'

  with pd.HDFStore(Mach15,'r') as d:
    dir = "shear_plots/"
    df = d.get("dat")
    dens_hot_Mach15 = np.array(df[df.columns[11]])
    vx_hot_Mach15 = np.array(df[df.columns[13]])
    vy_hot_Mach15 = np.array(df[df.columns[15]])
    vz_hot_Mach15 = np.array(df[df.columns[17]])
    vx_cl_Mach15 = np.array(df[df.columns[14]])
    vy_cl_Mach15 = np.array(df[df.columns[16]])
    vz_cl_Mach15 = np.array(df[df.columns[18]])
    size_Mach15 = np.array(df[df.columns[22]])
    rel_shear_Mach15 = np.sqrt((vx_hot_Mach15 - vx_cl_Mach15)**2 + (vy_hot_Mach15 - vy_cl_Mach15)**2 + (vz_hot_Mach15 - vz_cl_Mach15)**2.0)

    plt.hist(rel_shear_Mach15,6,density=True,facecolor='g',label=r"$\mathcal{M} = 1.5$")
    plt.xlabel("Shear Velocity")
    plt.ylabel("Probability")
    plt.grid(True)
    plt.legend()
    plt.savefig(dir+'plots/shear_pdf_'+str(i).zfill(3)+'.png')
    plt.close()


    plt.semilogy(rel_shear_Mach15, (np.array(size_Mach15))**(1./3.), 'go',label=r"$\mathcal{M} = 1.5$")
    plt.xlabel("Shear Velocity")
    plt.ylabel("Clump Radius")
    plt.ylim(0,100)
    plt.xlim(0,0.5)
    plt.legend()
    plt.savefig(dir+'plots/shear_size_'+str(i).zfill(3)+'.png')
    plt.close()

  with pd.HDFStore(Mach3,'r') as d:
    dir = "Mach3/"
    df = d.get("dat")
    dens_hot_Mach3 = np.array(df[df.columns[11]])
    vx_hot_Mach3 = np.array(df[df.columns[13]])
    vy_hot_Mach3 = np.array(df[df.columns[15]])
    vz_hot_Mach3 = np.array(df[df.columns[17]])
    vx_cl_Mach3 = np.array(df[df.columns[14]])
    vy_cl_Mach3 = np.array(df[df.columns[16]])
    vz_cl_Mach3 = np.array(df[df.columns[18]])
    size_Mach3 = np.array(df[df.columns[22]])
    rel_shear_Mach3 = np.sqrt((vx_hot_Mach3 - vx_cl_Mach3)**2 + (vy_hot_Mach3 - vy_cl_Mach3)**2 + (vz_hot_Mach3 - vz_cl_Mach3)**2.0)

    plt.hist(rel_shear_Mach3,6,density=True,facecolor='g',label=r"$\mathcal{M} = 3.0$")
    plt.xlabel("Shear Velocity")
    plt.ylabel("Probability")
    plt.grid(True)
    plt.legend()
    plt.savefig(dir+'plots/shear_pdf_'+str(i).zfill(3)+'.png')
    plt.close()


    plt.semilogy(rel_shear_Mach3, ((3./4.*np.pi)*np.array(size_Mach3))**(1./3.), 'go',label=r"$\mathcal{M} = 3.0$")
    plt.xlabel("Shear Velocity")
    plt.ylabel("Clump Radius")
    plt.ylim(0,100)
    plt.xlim(0,1.0)
    plt.legend()
    plt.savefig(dir+'plots/shear_size_'+str(i).zfill(3)+'.png')
    plt.close()

  with pd.HDFStore(Mach6,'r') as d:
    dir = "Mach6/"
    df = d.get("dat")
    dens_hot_Mach6 = np.array(df[df.columns[11]])
    vx_hot_Mach6 = np.array(df[df.columns[13]])
    vy_hot_Mach6 = np.array(df[df.columns[15]])
    vz_hot_Mach6 = np.array(df[df.columns[17]])
    vx_cl_Mach6 = np.array(df[df.columns[14]])
    vy_cl_Mach6 = np.array(df[df.columns[16]])
    vz_cl_Mach6 = np.array(df[df.columns[18]])
    size_Mach6 = np.array(df[df.columns[22]])
    rel_shear_Mach6 = np.sqrt((vx_hot_Mach6 - vx_cl_Mach6)**2 + (vy_hot_Mach6 - vy_cl_Mach6)**2 + (vz_hot_Mach6 - vz_cl_Mach6)**2.0)

    plt.hist(rel_shear_Mach6,6,density=True,facecolor='g',label=r"$\mathcal{M} = 6.0$")
    plt.xlabel("Shear Velocity")
    plt.ylabel("Probability")
    plt.grid(True)
    plt.legend()
    plt.savefig(dir+'plots/shear_pdf_'+str(i).zfill(3)+'.png')
    plt.close()


    plt.semilogy(rel_shear_Mach6, ((3./4.*np.pi)*np.array(size_Mach6))**(1./3.), 'go',label=r"$\mathcal{M} = 6.0$")
    plt.xlabel("Shear Velocity")
    plt.ylabel("Clump Radius")
    plt.ylim(0,100)
    plt.xlim(0,2.0)
    plt.legend()
    plt.savefig(dir+'plots/shear_size_'+str(i).zfill(3)+'.png')
    plt.close()

