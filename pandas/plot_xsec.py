from pandas import read_feather
import numpy as np
import matplotlib.pyplot as plt

import artie

NBINS = 50
EMIN  = 40
EMAX  = 70
MAX   = 1000000
TRES  = 125 # ns


# extract ideal cross-section measurement
df   = read_feather("ideal.feather")
tall = np.array(df.arrival_time/1000.0) # time 
eall = np.array(df.gen_energy*1000)     # energy 
edet = eall[(tall>0)]
hall,edges = np.histogram(eall, bins=NBINS,range=[EMIN, EMAX])
hdet,edges = np.histogram(edet, bins=NBINS,range=[EMIN, EMAX])
cbins = (edges[1:] + edges[:-1])/2.0

# obtain simulated data time stamps from MC:
tvac,evac = artie.simulate_data_from_mc("vacuum.feather")
targ,earg = artie.simulate_data_from_mc("argon.feather")


hvac,edges = np.histogram(evac, bins=NBINS,range=[EMIN, EMAX])
harg,edges = np.histogram(earg, bins=NBINS,range=[EMIN, EMAX])

# plot energy spectrum of detected vs produced neutrons:
plt.plot(cbins,hall,"b--",label="produced (generator)")
plt.plot(cbins,hdet,"r--",label="argon (ideal)")
plt.plot(cbins,hvac,"bo",label="vacuum (simulated)")
plt.plot(cbins,harg,"ro",label="argon (simulated)")
plt.xlabel("Kinetic Energy (keV)")
plt.ylabel("Neutrons")
plt.legend(frameon=False)
plt.savefig("ke.pdf")
plt.show()

# calculate transmission coefficient:
itr = np.zeros(cbins.size)
dtr = np.zeros(cbins.size)

mask = (hall > 0)
itr[mask] = hdet[mask] / hall[mask]
mask = (hvac > 0)
dtr[mask] = harg[mask] / hvac[mask]

plt.plot(cbins,itr,"k--",label="ideal")
plt.plot(cbins,dtr,"bo",label="simulation")
plt.legend(frameon=False)
plt.xlabel("Kinetic Energy (keV)")
plt.ylabel("Transmission Coefficient")
plt.legend(frameon=False)
plt.savefig("trans.pdf")
plt.show()

ixs = np.zeros(cbins.size)
dxs = np.zeros(cbins.size)
mask = (itr > 0)
ixs[mask] = - (1.0/4.2) * np.log(itr[mask])
mask = (dtr > 0)
dxs[mask] = - (1.0/4.2) * np.log(dtr[mask])
plt.plot(cbins,ixs,"k--",label="true")
plt.plot(cbins,dxs,"bo",label="simulation")
plt.xlabel("Kinetic Energy (keV)")
plt.ylabel("Cross Section (barns)")
plt.legend(frameon=False)
plt.savefig("xsec.pdf")
plt.show()









