import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt

def GaussianBeam(E0, w0, wavelength, x, y, z):
    """
    Define electric field of a gaussian beam, https://en.wikipedia.org/wiki/Gaussian_beam

    E0: electric field at the original
    w0: beam waist
    wavelength: laser wavelength
    x, y, z: coordinates, beam propagates along z-axis 
    """

    k = 2*np.pi/wavelength
    z_R = np.pi*(w0**2)/wavelength # Rayleigh range
    w_z = w0*np.sqrt(1+(z/z_R)**2)
    r = np.sqrt(x**2+y**2)
    inverse_Rz = z/(z**2+z_R**2) # inverse of R(z), radius of curvature
    phi_z = np.arctan(z/z_R) # Gouy phase

    return E0*w0/w_z*np.exp(-r**2/w_z**2)*np.exp(-1j*(k*z+k/2*inverse_Rz*r**2-phi_z))

def KirchoffIntegrand(E0, w0, wavelength, xp, yp, zp, x, y, z):
    k = 2*np.pi/wavelength
    R = np.sqrt((x-xp)**2+(y-yp)**2+(z-zp)**2)

    return 1/(1j*wavelength)*(np.exp(1j*k*R)/R)*((z-zp)/R)*GaussianBeam(E0, w0, wavelength, xp, yp, zp)



E0 = 1
w0 = 1e-3 # m
wavelength = 1e-6 # m
zp = 0
z = 0.1

x_list = np.linspace(-w0, w0, 2)
y_list = np.linspace(-w0, w0, 2)
intensity = np.zeros((len(x_list), len(y_list)))
for i, x in enumerate(x_list):
    for j, y in enumerate(y_list):
        intensity[i, j] = np.abs(dblquad(lambda xp, yp: KirchoffIntegrand(E0, w0, wavelength, xp, yp, zp, x, y, z), -2*w0, 2*w0, lambda x: -2*w0, lambda x: 2*w0)[0])**2

plt.imshow(intensity)
plt.show()