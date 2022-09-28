from pyexpat.model import XML_CQUANT_PLUS
import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt
import time
from multiprocessing import Pool

"""
To-do: try package xrt (https://xrt.readthedocs.io/index.html)
"""

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

def HollowGaussianBeam(E0, w0, wavelength, x, y, z, r):
    """
    Define electric field of a "hollow" gaussian beam

    E0: electric field at the original
    w0: beam waist
    wavelength: laser wavelength
    x, y, z: coordinates, beam propagates along z-axis 
    r: radius of hollow core of the beam
    """

    return GaussianBeam(E0, w0, wavelength, x, y, z)*(1-np.exp(-(x**2+y**2)/r**2))

def KirchhoffIntegrand(E0, w0, wavelength, xp, yp, zp, rp, x, y, z):
    """
    Define the integrand of Kirchhoff integral

    E0: electric field at the original
    w0: beam waist
    wavelength: laser wavelength
    xp, yp, zp: coordinates, where the electric field will be integrated, assume beam propagates along z-axis 
    rp: radius of hollow core of the incident beam
    x, y, z: coordinates, where the electric field is of interest, assume beam propagates along z-axis 
    """

    k = 2*np.pi/wavelength
    R = np.sqrt((x-xp)**2+(y-yp)**2+(z-zp)**2)

    return 1/(1j*wavelength)*(np.exp(1j*k*R)/R)*((z-zp)/R)*HollowGaussianBeam(E0, w0, wavelength, xp, yp, zp, rp)

def KirchhoffInteral(E0, w0, wavelength, zp, rp, x, y, z, integralrange):
    """
    Define Kirchhoff integral

    E0: electric field at the original
    w0: beam waist
    wavelength: laser wavelength
    zp: coordinate, where the electric field will be integrated, assume beam propagates along z-axis 
    rp: radius of hollow core of the incident beam
    x, y, z: coordinates, where the electric field is of interest, assume beam propagates along z-axis 
    integralrange: integral range, integrate from -integralrange/2 to integralrange/2

    It seems that the choice of z affects computation time, larger z leading to faster calculation.
    """

    realpart = dblquad(lambda xp, yp: np.real(KirchhoffIntegrand(E0, w0, wavelength, xp, yp, zp, rp, x, y, z)), 
                        -integralrange/2, integralrange/2, 
                        lambda x: -integralrange/2, lambda x: integralrange/2, 
                        epsrel = 1e-2)[0]

    imagpart = dblquad(lambda xp, yp: np.imag(KirchhoffIntegrand(E0, w0, wavelength, xp, yp, zp, rp, x, y, z)), 
                        -integralrange/2, integralrange/2, 
                        lambda x: -integralrange/2, lambda x: integralrange/2, 
                        epsrel = 1e-2)[0]

    r = realpart**2 + imagpart**2

    return r

def KirchhoffInteralMultiproc(E0, w0, wavelength, zp, rp, x_list, y, z, integralrange):
    """
    Define Kirchhoff integral

    E0: electric field at the original
    w0: beam waist
    wavelength: laser wavelength
    zp: coordinate, where the electric field will be integrated, assume beam propagates along z-axis 
    rp: radius of hollow core of the incident beam
    x_list, y, z: coordinates, where the electric field is of interest, assume beam propagates along z-axis 
    integralrange: integral range, integrate from -integralrange/2 to integralrange/2
    """

    E0_list = np.ones(len(x_list))*E0
    w0_list = np.ones(len(x_list))*w0
    wavelength_list = np.ones(len(x_list))*wavelength
    zp_list = np.ones(len(x_list))*zp
    rp_list = np.ones(len(x_list))*rp
    y_list = np.ones(len(x_list))*y
    z_list = np.ones(len(x_list))*z
    integralrange_list = np.ones(len(x_list))*integralrange

    param_list = list(zip(E0_list, w0_list, wavelength_list, zp_list, rp_list, x_list, y_list, z_list, integralrange_list))

    if __name__ == "__main__":
        with Pool(8) as p:
            intensity = np.array(p.starmap(KirchhoffInteral, param_list))

    return intensity


E0 = 1
w0 = 1e-3 # m, beam waist, = 2*sigma
wavelength = 1e-6 # m
zp = 0 # m
rp = w0/5
z = 0.2 # m
y = 0 # m
integralrange = 4*w0

x_list = np.linspace(0, 2*w0, 120)
numcycle = len(x_list)
counter = 1

if __name__ == "__main__":
    t0 = time.time()
    intensity = KirchhoffInteralMultiproc(E0, w0, wavelength, zp, rp, x_list, y, z, integralrange)
    print("Finishes in {:.2f} s".format(time.time()-t0))

    plt.plot(x_list/w0, intensity/(np.abs(E0)**2))
    plt.plot(x_list/w0, np.abs(HollowGaussianBeam(E0, w0, wavelength, x_list, y, z, rp)/E0)**2)
    plt.xlabel("x/$w_0$")
    plt.ylabel("Intensity/$E_0^2$")
    # plt.savefig("latest.jpg", dpi=600)
    plt.show()
