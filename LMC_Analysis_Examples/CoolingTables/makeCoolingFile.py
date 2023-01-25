import matplotlib
import numpy as np
from numpy import *
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate
import h5py

def keys(f):
        return [key for key in f.keys()]

filename = 'z_0.000.hdf5'
f = h5py.File(filename,'r+')

# List all groups
#print("Keys: %s" % keys(f))
a_group_key = keys(f)[12]  #12 is all metals
data = list(f[a_group_key])
#print(data[1])
#print(data)

temps = f['Total_Metals/Temperature_bins']
nH = f['Total_Metals/Hydrogen_density_bins']
net_cooling = f['Total_Metals/Net_cooling']
#print(net_cooling[:][:])

net_Cooling = net_cooling[:][:]
NH = nH[:]
Temps = temps[:]
print(net_Cooling)
##########################################
#net_heating = net_Cooling[net_Cooling < 0.0]
#net_Cooling = net_Cooling[net_Cooling > 0.0]

#separating heating from cooling
#Temps = Temps[net_Cooling > 0]
#temps_heat = Temps[net_Cooling < 0]
#NH = NH[net_Cooling > 0]
#nH_heat = NH[net_Cooling < 0]
#########################################


logtemps = zeros(len(Temps)-1)
lognH = zeros(len(NH)-1)
#logtemps_heat = zeros(len(temps_heat)-1)
#lognH_heat = zeros(len(nH_heat)-1)
lognet_cooling = [[0 for x in range(len(NH))] for x in range(len(Temps))]
#lognet_heating = [[0 for x in range(len(nH_heat))] for x in range(len(temps_heat))]

for i in np.arange(len(Temps)-1): 
    logtemps[i] = log10(Temps[i])
for i in np.arange(len(NH)-1): 
    lognH[i] = log10(NH[i])

for i in range(len(Temps)-1): 
    for j in range(len(NH)): 
       lognet_cooling[i][j] = log10(abs(net_Cooling[i][j]))
"""
for i in range(len(temps_heat)-1): 
    for j in range(len(nH_heat)): 
       lognet_heating[i][j] = log10(abs(net_heating[i][j]))
"""
dens,tmp = np.meshgrid(lognH,logtemps)
#dens_heat,tmp_heat = np.meshgrid(lognH_heat,logtemps_heat)

#lognet_cooling = lognet_cooling.reshape(dens.shape)

#xi = np.linspace(min(logtemps), max(logtemps))
#yi = np.linspace(min(lognH), max(lognH))

#print(temps[20])
#print(nH[:])
#lognet_cooling = griddata(logtemps,lognH,lognet_cooling,xi,yi)

cooling_nH1 = zeros(len(Temps)-1)
cooling_nH2 = zeros(len(Temps)-1)
cooling_nH3 = zeros(len(Temps)-1)
cooling_nH1_nolog = zeros(len(Temps)-1)
cooling_nH2_nolog = zeros(len(Temps)-1)
cooling_nH3_nolog = zeros(len(Temps)-1)

#heating_nH1 = zeros(len(temps_heat)-1)
#heating_nH2 = zeros(len(temps_heat)-1)

print(nH[len(nH)-20])
print(nH[len(nH)-50])
#print(nH_heat[len(nH_heat)-20])
#print(nH_heat[len(nH_heat)-50])


for i in np.arange(len(temps[:])-1):
   # print(net_cooling[i][len(NH)-50])
    cooling_nH1[i] = lognet_cooling[i][len(NH)-1]
    cooling_nH2[i] = lognet_cooling[i][len(NH)-20]
    cooling_nH3[i] = lognet_cooling[i][len(NH)-50]
    cooling_nH1_nolog[i] = -net_cooling[i][len(NH)-1]
    cooling_nH2_nolog[i] = -net_cooling[i][len(NH)-20]
    cooling_nH3_nolog[i] = -net_cooling[i][len(NH)-50]

"""
plt.plot(logtemps,cooling_nH1_nolog,'b',label='nH = 1')
plt.plot(logtemps,cooling_nH2_nolog,'b--',label='nH = 1e-2')
plt.plot(logtemps,cooling_nH3_nolog,'b-.',label='nH = 1e-5')
plt.xlabel("log(T) (K)")
plt.ylim(-1e-21,0)
plt.ylabel(r"$\Lambda$")
plt.legend()
#plt.savefig("Cooling.pdf")
plt.close()
"""

"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(dens,tmp,lognet_cooling[:][:])
#ax.set_xscale("log", nonposx='clip')
#ax.set_zscale("log")
#plt.xlim(-5,0)
ax.set_zlim(-24,-21)
ax.set_xlim(-5,0)
ax.set_ylim(3,8)
plt.xlabel('Hydrogen Number Density')
plt.ylabel('Temperature (K)')
plt.title('|Normalized Lambda|')
plt.show()
"""

"""
plt.contourf(dens,tmp,net_cooling[:],norm=matplotlib.colors.LogNorm())
plt.xlabel('Hydrogen Number Density')
plt.ylabel('Temperature (K)')
plt.title('Normalized Lambda')
plt.colorbar()
plt.show()
"""
#print(temps[:])

data = loadtxt('cooling_abund03.txt')
logT,kt,totlam,lam = data[:,0],data[:,1],data[:,2],data[:,3]

loglam = zeros(len(logT))
totloglam = zeros(len(logT))
for i in np.arange(len(logT)):
    loglam[i] = log10(lam[i])
    totloglam[i] = log10(totlam[i])

f = open('analyticCooling.txt', 'r')
logTSchure = []
analytic = []
lamSchure = []

while True:
    s = f.readline()
    if not s: break
    logTSchure.append(float(s[2:7]))
    analytic.append(float(s[7:20]))
    lamSchure.append(float(s[20:34]))

trueAnalytic = zeros(len(logTSchure))
T = zeros(len(logTSchure))
theta = zeros(len(logTSchure))
lamHensler = zeros(len(logTSchure))
theta03 = zeros(len(logTSchure))
thetatown = zeros(len(logTSchure))
z3fit = zeros(len(logTSchure))
townsleyFit = zeros(len(logTSchure))

T0 = 2e5
T02 = 1.5e5
T03 = 3.0e5

#######################
# Hensler 1987 fit
z = 0.0 #metallicity compared to solar
A5 = (2.5+7.0*(z**0.5))/(5.0-log10(1.48e11*(z**1.1)+1.0e6))
A4 = 30.0*z
lambda05 = 10.0**(-22.0-(5.0*A5)+7.0*(z**0.5))
lambda04 = 10.0**(-22.0-(4.0*A4))
print(lambda05)
print(lambda04)
############################

###########################
# Photoionization Heating
#
logheat = zeros(len(logTSchure))
T_background = 44000
psi = 1.52
phi2 = 0.7
xi2 = 0.38
KB = 1.3806e-16
for i in np.arange(len(logTSchure)-1):
        T[i] = 10**logTSchure[i]
        logheat[i] = log10(abs(2.07*10**(-11)*(psi*KB*T_background*phi2 - KB*T[i]*xi2)/sqrt(T[i])))


for i in np.arange(len(logTSchure)-1):
    T[i] = 10**logTSchure[i]
    theta[i] = 0.4*log10(T[i]/T0) - 3 + 6.2/(exp(1.5*log10(T[i]/T0)+0.08) + exp(-(log10(T[i]/T0) + 0.08)))
    trueAnalytic[i] = log10((1.5e-21)*10**theta[i])
    if (T[i] >= 10**5):
        lamHensler[i] = log10(lambda05*(T[i]**A5))
    else:   
        lamHensler[i] = log10(lambda04*(T[i]**A4))  
    theta03[i] = 0.4*log10(T[i]/T02) - 3 + 4.6/(exp(2.2*log10(T[i]/T02)+0.08) + exp(-1.0*(log10(T[i]/T02) + 0.08)))
    thetatown[i] = 0.4*log10(T[i]/T02) - 3 + 5.2/(exp(2.9*log10(T[i]/T03)+0.08) + exp(-1.0*(log10(T[i]/T03) + 0.08)))
    z3fit[i] = log10((1.5e-21)*10**theta03[i])
    townsleyFit[i] = log10((1.5e-21)*10**thetatown[i])

plt.plot(logtemps,cooling_nH1,'r',label='Cloudy, nH = 1')
plt.plot(logtemps,cooling_nH2,'r--',label='Cloudy, nH = 1e-2')
plt.plot(logtemps,cooling_nH3,'r-.',label='Cloudy, nH = 1e-5')
plt.plot(logTSchure,lamSchure,'k',label="Z=1.0 (Schure+ 2009)")
#plt.plot(logTSchure,trueAnalytic,'g--',label = "approximation to Z=1.0 (BuZD2016)")
#plt.plot(logTSchure,z3fit,'r--',label = "fit to Z=0.3 (Me)")
#plt.plot(logTSchure,townsleyFit,'r--',label = "approximation to 30 Dor (this work))",linewidth=2)
#plt.plot(logTSchure,lamHensler,'g',label = "Z=0.3 approx (Hensler 1987)")
#plt.plot(logT,totloglam,'b--',label= "Z=0.3; 0.1-50 keV (XSPEC)")
#plt.plot(logT,loglam,'b',label= "Z=0.3; 0.5-8 keV (XSPEC)")
plt.plot(logTSchure,logheat,'b',label= "Spitzer")
plt.xlabel('log Temperature (K)',fontsize=16)
plt.ylabel(r'log($\Lambda(T)$) (ergs/s $cm^{-3}$)',fontsize=16)
plt.xlim(3.5,8)
plt.ylim(-24.5,-19.5)
plt.legend(loc=1,prop={'size':11.5})
plt.title('Cooling Curves')
#plt.show()
plt.savefig('photoionization_cloudy_lambda.pdf')
plt.close()
"""
#logT = np.arange(3.8,8.2,0.04)
#kt = zeros(len(logT))
#for i in range(len(logT)):
#    kt[i] = (8.6e-8)*10**logT[i]

#np.savetxt('cooling_30DorFit.txt',np.c_[logTSchure,],fmt='%.3e')
f = open('cooling_30Dor.txt','w')
for i in range(len(logT)):
    f.write('{:12f}'.format(logT[i]))
    f.write('{:12f}'.format(totloglam[i]))
    f.write('{:12f}'.format(loglam[i]))
    f.write("\n")

f.close()

"""
