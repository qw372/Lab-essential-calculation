import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# only applies at low intensity, ie, I << Isat
f87 = (384.230484468-0.002563006+0.000193741)*1e12  # Hz, Rb87 D2 F=2->F'=3
f85 = (384.230406373-0.001264889+0.000100205)*1e12  # Hz, Rb85 D2 F=3->F'=4
isotope_shift = f85-f87  # Hz

h = 6.62607015e-34 # Joule x second
c = 299792458 # m/s
Gamma = 6.0666e6*2*np.pi # Hz, Rb87 D2 linen natural linewidth
wavevec87 = 2*np.pi*f87/c # 1/meter, laser wavevector
wavevec85 = 2*np.pi*f85/c # 1/meter, laser wavevector
mass87 = 0.001/(6.022e23)*87 # kg, Rb87 mass
mass85 = 0.001/(6.022e23)*85 # kg, Rb87 mass
cross_section0 = (780.24e-9)**2/2/np.pi # meter^2, Rb87(85), D2, lambad^2/2/Pi
Kb = 1.38e-23 # m^2*kg/s^2/Kelvin

Rb87e_F = np.array([3, 2, 1, 0]) # Rb87 excited state F quantum number
Rb87g_F = np.array([2, 1])
Rb87e_split = np.array([0, -266.65e6, (-266.65-156.947)*1e6, (-266.65-156.947-72.2180)*1e6])
Rb85e_F = np.array([4, 3, 2, 1])
Rb85g_F = np.array([3, 2])
Rb85e_split = np.array([0, -120.64e6, (-120.64-63.401)*1e6, (-120.64-63.401-29.372)*1e6])

eta87 = 0.2783 # natural abundance of Rb87
eta85 = 0.7217 # natural abundance of Rb85

def integrand87(v, T, detuning):
    return np.exp(-mass87*v**2/2/Kb/T)/(1+4*(2*np.pi*detuning-wavevec87*v)**2/Gamma**2)
def integrand85(v, T, detuning):
    return np.exp(-mass85*v**2/2/Kb/T)/(1+4*(2*np.pi*(detuning-isotope_shift)-wavevec85*v)**2/Gamma**2)
# detuning here refers to laser detuning w.r.t. Rb87 D2 F=2->F'=3 transition

def crossRb87(T, detuning):
    cs = np.sum([cross_section0*(2*Rb87e_F[i]+1)/(2*Rb87g_F[0]+1)*np.sqrt(mass87/2/np.pi/Kb/T)
                                                      *integrate.quad(integrand87, -2000, 2000, args=(T, detuning-Rb87e_split[i]))[0]
                                                      for i in np.arange(len(Rb87e_F)-1)])
    return cs
# returns Rb87 D2 F=2->F'=1/2/3 transition cross section
# presume that Rb population uniformly distributes among ground hyperfine states

def crossRb85(T, detuning):
    cs = np.sum([cross_section0*(2*Rb85e_F[i]+1)/(2*Rb85g_F[0]+1)*np.sqrt(mass85/2/np.pi/Kb/T)
                                                      *integrate.quad(integrand85, -2000, 2000, args=(T, detuning-Rb85e_split[i]))[0]
                                                      for i in np.arange(len(Rb85e_F)-1)])
    return cs 
# returns Rb85 D2 F=3->F'=2/3/4 transition cross section
# presume that Rb population uniformly distributes among ground hyperfine states

def absorption(T, vaporpres, L, detuning):
    N = vaporpres/Kb/T # Rb density based on ideal gas law
    NRb87 = N*eta87*(Rb87g_F[0]*2+1)/np.sum(Rb87g_F*2+1) # F=1 ground state Rb87 density
    NRb85 = N*eta85*(Rb85g_F[0]*2+1)/np.sum(Rb85g_F*2+1) # F=2 ground state Rb85 density
    return np.exp(-L*(crossRb87(T, detuning)*NRb87+crossRb85(T, detuning)*NRb85))