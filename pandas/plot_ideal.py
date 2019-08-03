from pandas import read_feather
import numpy as np
import matplotlib.pyplot as plt

NBINS = 200
EMIN  = 40
EMAX  = 70
MAX   = 1000000
TRES  = 125 # ns

df = read_feather("ideal.feather")
tall = np.array(df.arrival_time/1000.0) # time 
eall = np.array(df.gen_energy*1000)     # energy 
epss = eall[(tall>0)]

hall,edges = np.histogram(eall, bins=NBINS,range=[EMIN, EMAX])
hpss,edges = np.histogram(epss, bins=NBINS,range=[EMIN, EMAX])
cbins = (edges[1:] + edges[:-1])/2.0

# plot energy spectrum of detected vs produced neutrons:
plt.plot(cbins,hall,"b--",label="produced")
plt.plot(cbins,hpss,"r--",label="detected")
plt.xlabel("Kinetic Energy (keV)")
plt.ylabel("Neutrons")
plt.legend(frameon=False)
plt.savefig("itke.pdf")
plt.show()

# calculate transmission coefficient:
tr = np.zeros(cbins.size)
mask = (hall > 0)
tr[mask] = hpss[mask] / hall[mask]

plt.plot(cbins,tr,"k--",label="true")
plt.savefig("itrans.pdf")
plt.show()

xr = np.zeros(tr.size)
mask = (tr > 0)
xr[mask] = - (1.0/4.2) * np.log(tr[mask])
plt.plot(cbins,xr,"k--",label="true")
plt.show()









