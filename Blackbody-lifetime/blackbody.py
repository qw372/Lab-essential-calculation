import numpy as np
from sympy.physics.wigner import *

def BBRSpectrum(nu, T):
    """
    Blackbody radiationi spectrum
    nu: frequency in 1/cm
    T: temperature in Kelvin
    """

    h = 6.62607015e-34 # SI units, Planck's constant
    c = 299792458 # SI units, speed of light
    kB = 1.380649e-23 # SI units, Boltzmann constant

    nu *= (100*c) # convert nu to Hz

    return 8*np.pi*h*nu**3/c**3/(np.exp(h*nu/kB/T)-1)


