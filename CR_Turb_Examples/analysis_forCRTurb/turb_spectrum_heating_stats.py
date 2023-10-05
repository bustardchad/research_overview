import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units.yt_array import YTQuantity
from yt.fields.api import ValidateParameter
import scipy.stats as stats
pUnit = YTQuantity(1, 'cm**2/s**2')
PEUnit = YTQuantity(1, 'cm/s**2')

Myr = 1.

vmax = 5.0

#Streaming energy loss as a function of fields output from Athena++ simulation
def heating(field,data):
          vb = (data['velocity_x']*data['magnetic_field_x'] + data['velocity_y']*data['magnetic_field_y'] + data['velocity_z']*data['magnetic_field_z'])*YTQuantity(1,"s/cm")/data['magnetic_field_magnitude']
          sigma1 = (data['alfven_speed']*YTQuantity(1,"s/cm"))*(1./(1./data['Sigma_adv1'] + 1./data['Sigma_diff1']))
          s1 = data['Fc1'] - (4./(3.))*data['Ec']*vb/vmax
          return np.abs(sigma1*s1)*8.0/(3e-5)


yt.add_field(("gas","heating"), function=heating,display_name="Streaming Energy Loss Rate",units="")


def create_spectrum_binned_statistic(ds,i):

    # a FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this.

    max_level = ds.index.max_level

    low = ds.domain_left_edge
    dims = ds.domain_dimensions*int(1.0)
    nx, ny, nz = dims
    
    nindex_rho = 1./2.
    Kk = np.zeros( (nx//2+1, ny//2+1, nz//2+1))
   
    Kk = fft_comp(ds, ("gas","heating"), max_level, low, dims)
    # wavenumbers
    L = (ds.domain_right_edge - ds.domain_left_edge).d
    kx = np.fft.fftfreq(nx)*nx/L[0]
    ky = np.fft.fftfreq(ny)*ny/L[1]
    kz = np.fft.fftfreq(nz)*nz/L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0/L)
    kmax = np.min(0.5*dims/L)

    kbins = np.arange(0.5, nx//2+1, 1.)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, kx, kx)
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)
    knrm = k.flatten()
    E_spectrum, _, _ = stats.binned_statistic(knrm, Kk.flatten(),
                                      statistic = "mean",
                                      bins = kbins)
    E_spectrum *= 4.*np.pi/3. * (kbins[1:]**3 - kbins[:-1]**3)
    
    k = 0.5 * (kbins[1:] + kbins[:-1])


    
    index = np.argmax(E_spectrum)
    kmax = k[index]
    print("Wavelength with highest energy density (in kpc): ")
    print(1.0/(kmax))
    Emax = E_spectrum[index]
    print("Emax: ")
    print(Emax)

    print("Energy spectrum: ")
    print(E_spectrum)

    return k, E_spectrum

def fft_comp(ds, iu, level, low, delta ):

    cube = ds.covering_grid(level, left_edge=low,
                            dims=delta,
                            fields=[iu])

    u = cube[iu].d
    ru = np.fft.fftn(u)

    return np.abs(ru)**2  # rho v^2 




ts = yt.DatasetSeries("../cr.out1.0007*")
i = 0
EkinVector = []
EkinrhoVector = []
EthermVector = []
EthermrhoVector = []
magEVector = []
timeVector = []
ratioVector = []
ratiorhoVector = []
totalEVector = []
gravPEVector = []
log10EkinrhoVector = []
teddyvec = []
ekvec = np.zeros(128)

for ds in ts:
    dd = ds.all_data()
    teddyvec.append(ds.current_time.v*3.155e13/(3.0856e21*0.666667/5e6))
    timeVector.append(ds.current_time.v/Myr)
    kvec,ekvecout = create_spectrum_binned_statistic(ds,i)
    ekvec = np.vstack([ekvec, ekvecout]) # does the eigenmode analysis and spits out KE vs k
    i = i+1

a = ekvec[1:9,:]  # last few rows (times) of spectrum array

avg = np.mean(a,axis = 0)
mina = np.min(a,axis = 0)
maxa = np.max(a,axis = 0)


k = kvec[0:len(avg)]


print("pcr/pg = 1, beta = 1, res = 256: ")
print(avg)
print(mina)
print(maxa)
print("k: ")
print(k)

plt.loglog((k), avg*(k**2), 'bo-',label=r"$P_{CR}/P_{g} \sim 1, \beta \sim 1$")
plt.fill_between((k), mina*(k**2), maxa*(k**2), facecolor='blue', alpha=0.5)
plt.ylim(1E-5,4E-3)
#plt.xlim(3E-2,2E0)
plt.xlim(0.5,50)
#plt.xlabel(r"$\lambda / L$",fontsize=18)
plt.xlabel(r"$k$",fontsize=18)
plt.ylabel(r"FFT of CR Heating Rate",fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig("heating_spectrum_streaming_res256_pcpg1_varybeta.pdf")
plt.close()
