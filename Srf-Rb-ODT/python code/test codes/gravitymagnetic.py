import numpy as np
import matplotlib.pyplot as plt

def trap_gravity(x, optical_trap_depth, sigma, m, Bfieldgrad, magneton):
    # x is position, in um
    # optical_trap_depth in uK
    # sigma is sigma of ODT beam gaussian distribution, in um
    # m is mass in atomic mass unit
    # Bfieldgrad in G/cm
    # magneton in MHz/G

    kB = 1.380649e-23 # SI units
    x *= 1.e-6 # convert to m
    sigma *= 1.e-6 # convert to m
    m *= 1.66053906660e-27 # kg
    g = 9.8 # kg*m/s^2
    h = 6.62607015e-34 # SI unit

    r = -(optical_trap_depth/1.e6)*kB*np.exp(-x**2/2/sigma**2) + m*g*x + h*(Bfieldgrad*1e2*1e6)*magneton*x

    return r/kB*1e6

def trap_gravity_derivative(x, optical_trap_depth, sigma, m, Bfieldgrad, magneton):
    # x is position, in um
    # optical_trap_depth in uK
    # sigma is sigma of ODT gaussian distribution, in um
    # m is mass in atomic mass unit

    kB = 1.380649e-23 # SI units
    x = x*1.e-6 # convert to m, x is an array
    sigma *= 1.e-6 # convert to m
    m *= 1.66053906660e-27 # kg
    g = 9.8 # kg*m/s^2

    r = (optical_trap_depth/1.e6)*kB*np.exp(-x[0]**2/2/sigma**2)*(x[0]/sigma**2) + m*g + Bfieldgrad*magneton

    return r/kB

m_Rb = 86.9 # dalton
dalton = 1.66053906660e-27 # kg
g = 9.8 # m/s^2
h = 6.62607015e-34 # SI unit
gravity_Rb = m_Rb*dalton*g/h/1e6/1e2 # MHz/cm
bohr_magneton = 1.399624624 # MHz/G
g_factor_Rb = 0.5
mF_Rb = 1
magneton_Rb = g_factor_Rb*bohr_magneton*mF_Rb # MHz/G

gravity_Bgradient = gravity_Rb/magneton_Rb # B-field gradient required to balance gravity


optical_trap_depth = 22 # uK
sigma = 20 # um
Bfieldgrad = 140 # G/cm
x = np.linspace(-100, 100, 1000)

y = [trap_gravity(i, optical_trap_depth, sigma, m_Rb, Bfieldgrad, magneton_Rb) for i in x]
plt.plot(x, y)
plt.grid()
plt.show()