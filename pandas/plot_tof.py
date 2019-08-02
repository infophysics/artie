from pandas import read_feather
import numpy as np
import matplotlib.pyplot as plt

def energy_from_tof(tof):
    MN = 939*1000.0  # keV 
    C  = 3E2         # m / mus
    L = 70        # m
    beta = (L / tof) / C
    gamma = (1.0 / (1 - beta**2))**0.5
    return (gamma - 1.0)*MN    

df = read_feather("ideal.feather")
tall = np.array(df.arrival_time/1000.0) # time 
eall = np.array(df.gen_energy*1000)     # energy 
tpss = tall[(tall>0)]
epss = eall[(tall>0)]
etof = energy_from_tof(tpss)

tfin  = np.linspace(18,48,200)
efin = energy_from_tof(tfin)

plt.plot(epss, tpss, "b.", label="GEANT")
plt.plot(efin, tfin, "r--", label="Einstein")
plt.xlim(40,70)
plt.xlabel("Kinetic Energy (keV)")
plt.ylabel("Time-of-Flight (mu s)")
plt.legend(frameon=False,loc=1)
plt.savefig("itof.pdf")
plt.show()

h,edges = np.histogram(etof-epss , bins=200,range=[-10.0,10.0])
cbins = (edges[:-1] + edges[1:])/2.0
plt.plot(cbins,h,"b-")
plt.xlabel("Kinetic Energy (measured - actual) (keV)")
plt.savefig("idt.pdf")
plt.show()




