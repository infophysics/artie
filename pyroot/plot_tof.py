import ROOT
import numpy as np
import matplotlib.pyplot as plt

def energy_from_tof(tof):
    MN = 939*1000.0  # keV 
    C  = 3E2         # m / mus
    L = 100        # m
    beta = (L / tof) / C
    gamma = (1.0 / (1 - beta**2))**0.5
    return (gamma - 1.0)*MN    
    
fv = ROOT.TFile("vacuum.root")
tv = fv.Get("artie")
#tv.Print()

fa = ROOT.TFile("argon.root")
ta = fa.Get("artie")
ta.Print()

tar = np.array([])
ear = np.array([])

for evt in ta:
    t = evt.arrival_time/1000.0  # mus
    e = evt.gen_energy*1000.0    # keV
    if (t > 0):
        #print(t, e)
        tar = np.append(tar, t)
        ear = np.append(ear, e)



tfin  = np.linspace(18,48,200)
efin = energy_from_tof(tfin)

plt.plot(ear,  tar, "b.", label="GEANT")
plt.plot(efin, tfin, "r--", label="Einstein")
plt.xlim(40,70)
plt.xlabel("Kinetic Energy (keV)")
plt.ylabel("Time-of-Flight (mu s)")
plt.legend(frameon=False)
plt.savefig("tof_kvt_true.pdf")
plt.show()

h,edges = np.histogram(energy_from_tof(tar) - ear, bins=200,range=[-10.0,10.0])
cbins = (edges[1:] + edges[:-1])/2.0

plt.plot(cbins,h,"b-")
plt.xlabel("Kinetic Energy (measured - actual) (keV)")
plt.savefig("tof_dt_true.pdf")
plt.show()

tar = tar + np.random.normal(size=tar.size)*0.120

plt.plot(ear,  tar, "b.", label="GEANT")
plt.plot(efin, tfin, "r--", label="Einstein")
plt.xlim(40,70)
plt.xlabel("Kinetic Energy (keV)")
plt.ylabel("Time-of-Flight (mu s)")
plt.legend(frameon=False)
plt.savefig("tof_kvt_meas.pdf")
plt.show()

h,edges = np.histogram(energy_from_tof(tar) - ear, bins=200,range=[-10.0,10.0])
cbins = (edges[1:] + edges[:-1])/2.0

plt.plot(cbins,h,"b-")
plt.xlabel("Kinetic Energy (measured - actual) (keV)")
plt.savefig("tof_dt_meas.pdf")
plt.show()





