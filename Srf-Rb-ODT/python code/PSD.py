from matplotlib.pyplot import psd
import numpy as np

T = 100e-9 # K
n = 2e9*1e6 # density, 1/m^3
n *= 1700/9 # crossed-beam ODT
n *= 20 # gain from better slowing scheme

hbar = 1.05457182e-34 # SI units
m = 106.62*1.66053906660e-27 # kg
kB = 1.380649e-23 # J/K

thermal_wavelength = hbar*np.sqrt(2*np.pi/m/kB/T)
PSD = (thermal_wavelength**3)*n

print(PSD)