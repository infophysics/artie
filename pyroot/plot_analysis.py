import ROOT
import numpy as np
import matplotlib.pyplot as plt

def energy_from_tof(tof):
    MN = 939*1000.0  # keV 
    C  = 3E2         # m / mus
    L = 71.03        # m
    beta = (L / tof) / C
    gamma = (1.0 / (1 - beta**2))**0.5
    return (gamma - 1.0)*MN    
    
fv = ROOT.TFile("vacuum.root")
tv = fv.Get("artie")
#tv.Print()

fa = ROOT.TFile("argon.root")
ta = fa.Get("artie")
ta.Print()

tar = np.array([]) # time argon
ear = np.array([]) # energy argon
tva = np.array([]) # time vacuum
eva = np.array([]) # energy vacuum

for evt in ta:
    t = evt.arrival_time/1000.0  # mus
    e = evt.gen_energy*1000.0    # keV
    if (t > 0):
        #print(t, e)
        tar = np.append(tar, t)
        ear = np.append(ear, e)

for evt in tv:
    t = evt.arrival_time/1000.0  # mus
    e = evt.gen_energy*1000.0    # keV
    if (t > 0):
        #print(t, e)
        tva = np.append(tva, t)
        eva = np.append(eva, e)

tar = tar + np.random.normal(size=tar.size)*0.120
tva = tva + np.random.normal(size=tva.size)*0.120

ema  = energy_from_tof(tar) # energy measured argon
emv  = energy_from_tof(tva) # energy measured vacuum

NBINS = 200

hear,edges = np.histogram(ear, bins=NBINS,range=[40.0,70.0])
heva,edges = np.histogram(eva, bins=NBINS,range=[40.0,70.0])
cbins = (edges[1:] + edges[:-1])/2.0
plt.plot(cbins,hear,"b-",label="argon")
plt.plot(cbins,heva,"r--",label="vacuum")
plt.xlabel("True Kinetic Energy (keV)")
plt.ylabel("Neutrons Detected")
plt.legend(frameon=False)
plt.savefig("tke.pdf")
plt.show()


hema,edges = np.histogram(ema, bins=NBINS,range=[40.0,70.0])
hemv,edges = np.histogram(emv, bins=NBINS,range=[40.0,70.0])
cbins = (edges[1:] + edges[:-1])/2.0
plt.plot(cbins,hema,"b-",label="argon")
plt.plot(cbins,hemv,"r--",label="vacuum")
plt.xlabel("Measured Kinetic Energy (keV)")
plt.ylabel("Neutrons Detected")
plt.legend(frameon=False)
plt.savefig("mke.pdf")
plt.show()

tr = np.zeros(cbins.size)
mask = (heva > 0)
tr[mask] = hear[mask] / heva[mask]

tm = np.zeros(cbins.size)
fm = np.zeros(cbins.size)
mask = (hemv > 0)
tm[mask] = hema[mask] / hemv[mask]
fm[mask] = np.sqrt(1.0/hemv[mask] + 1.0/hema[mask])
plt.plot(cbins,tr,"r--",label="true")
plt.errorbar(cbins,tm,yerr=tm*fm,label="measured",fmt="bo")
plt.xlabel("Kinetic Energy (keV)")
plt.ylabel("Transmission Coefficient")
plt.legend(frameon=False)
plt.savefig("trans.pdf")
plt.show()

xr = np.zeros(tr.size)
mask = (tr > 0)
xr[mask] = - (1.0/4.2) * np.log(tr[mask])

xm = np.zeros(tm.size)
sm = np.zeros(tm.size)
mask = (tm > 0)
xm[mask] = - (1.0/4.2) * np.log(tm[mask])
sm[mask] = 10*(1.0/4.2) * sm[mask] / tm[mask]

plt.plot(cbins,xr,"r--",label="true")
plt.errorbar(cbins,xm,yerr=xm*fm,fmt="bo",label="measured")
plt.xlabel("Kinetic Energy (keV)")
plt.ylabel("Cross Section")
plt.yscale("log")
plt.savefig("xsec.pdf")
plt.show()










