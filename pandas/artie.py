# useful routines for Artie analysis

from pandas import read_feather
import numpy as np

def energy_from_tof(tof):
    MN = 939.5654133*1000.0  # keV 
    C  = 2.99792458*100.0    # m / mus
    L = 70        # m
    beta = (L / tof) / C
    gamma = (1.0 / (1 - beta**2))**0.5
    return (gamma - 1.0)*MN    

def smear_time(t):
    RES = 0.125  # time resolution in mus
    t = t + np.random.uniform(-RES/2.0, RES/2.0, size=t.size) 
    return t

def simulate_data_from_mc(filename):
    df = read_feather(filename)
    t = np.array(df.arrival_time/1000.0)
    t = t[(t > 0)]
    t = smear_time(t)
    e = energy_from_tof(t)
    return t,e
