from pandas import read_feather
import numpy as np
import matplotlib.pyplot as plt
import artie

df = read_feather("ideal.feather")
tall = np.array(df.arrival_time/1000.0) # time 
eall = np.array(df.gen_energy*1000)     # energy 
tpss = tall[(tall>0)]
epss = eall[(tall>0)]
etof = artie.energy_from_tof(tpss)

tfin  = np.linspace(15,30,150)
efin = artie.energy_from_tof(tfin)

plt.plot(epss, tpss, "b.", label="GEANT")
plt.plot(efin, tfin, "r--", label="Einstein")
plt.xlim(40,70)
plt.xlabel("Kinetic Energy (keV)")
plt.ylabel("Time-of-Flight (mu s)")
plt.legend(frameon=False,loc=1)
plt.savefig("itof.pdf")
plt.show()

delta = etof - epss
print("mean time difference:  ", np.mean(delta))
print("variance:              ", np.var(delta))
h,edges = np.histogram(etof-epss , bins=201,range=[-10.0,10.0])
cbins = (edges[:-1] + edges[1:])/2.0
plt.plot(cbins,h,"b-")
plt.xlabel("Kinetic Energy (measured - actual) (keV)")
plt.savefig("idt.pdf")
plt.show()




