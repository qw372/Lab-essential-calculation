import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from multiprocessing import Pool

wavelength = 20 # mm, microwave wavelength
k = 2*np.pi/wavelength
Omega0 = 2*np.pi*35 # kHz, peak Rabi frequency
sigma_x = 1.1 # mm, ODT cloud size
sigma_y = sigma_x
sigma_z = sigma_x

def damped_oscillation(t, Omega, tau):
    return 0.5 + 0.5 * np.cos(Omega*t) * np.exp(-t/tau)

def integrated_signal_1D(t_list, alpha, x0):

    rabi_freq = lambda x: Omega0 * (1 + alpha*np.cos(2*k*x)) / ((1 + alpha*np.cos(2*k*x0)))
    Rabi_oscillation = lambda t, x: 0.5 + 0.5 * np.cos(rabi_freq(x) * t)
    density_dist = lambda x: np.exp(-0.5*((x-x0)/sigma_x)**2) / (sigma_x*(2*np.pi)**(1/2)) # gaussian

    signal = np.empty(len(t_list))
    for i, t in enumerate(t_list):
        integrand = lambda x: Rabi_oscillation(t, x) * density_dist(x)
        signal[i] = integrate.quad(integrand, x0-sigma_x*5, x0+sigma_x*5, epsabs=1e-3, epsrel=1e-3)[0]

    return signal

def eval_signal_1D(t_list, alpha_list, x0_list):
    signal = np.empty((len(alpha_list), len(x0_list), len(t_list)))
    for i, alpha in enumerate(alpha_list):
        for j, x0 in enumerate(x0_list):
            signal[i, j] = integrated_signal_1D(t_list=t_list, alpha=alpha, x0=x0)

    return signal

def integrated_signal_3D(t_list, alpha, x0, y0, z0):

    rabi_freq = lambda x, y, z: Omega0 * (1 + alpha*np.cos(2*k*x)) / ((1 + alpha*np.cos(2*k*x0))) * (1 + alpha*np.cos(2*k*y)) / ((1 + alpha*np.cos(2*k*y0))) * (1 + alpha*np.cos(2*k*z)) / ((1 + alpha*np.cos(2*k*z0)))
    Rabi_oscillation = lambda t, x, y, z: 0.5 + 0.5 * np.cos(rabi_freq(x, y, z) * t)
    density_dist = lambda x, y, z: np.exp(-0.5*((x-x0)/sigma_x)**2 - 0.5*((y-y0)/sigma_y)**2 - 0.5*((z-z0)/sigma_z)**2) / (sigma_x*sigma_y*sigma_z*(2*np.pi)**(3/2)) # gaussian

    signal = np.empty(len(t_list))
    for i, t in enumerate(t_list):
        integrand = lambda x, y, z: Rabi_oscillation(t, x, y, z) * density_dist(x, y, z)
        signal[i] = integrate.tplquad(integrand, z0-sigma_z*3, z0+sigma_z*3, y0-sigma_y*3, y0+sigma_y*3, x0-sigma_x*3, x0+sigma_x*3, epsabs=1e-2, epsrel=1e-2)[0]

    return signal

def eval_signal_3D(t_list, alpha_list, x0_list):
    alpha_list_repeated = np.repeat(alpha_list, len(x0_list))
    x0_list_repeated = np.tile(x0_list, len(alpha_list))
    t_list_repeated = np.tile(t_list, len(alpha_list)*len(x0_list)).reshape(len(alpha_list)*len(x0_list), len(t_list))
    param_list = list(zip(t_list_repeated, alpha_list_repeated, x0_list_repeated, x0_list_repeated, x0_list_repeated))

    if __name__ == "__main__":
        with Pool(12) as p:
            signal = np.array(p.starmap(integrated_signal_3D, param_list))

            signal = np.reshape(signal, (len(alpha_list), len(x0_list), len(t_list)))

    return signal

if __name__ == "__main__":
    x0_list = [0, wavelength/16, wavelength/8, wavelength*3/16, wavelength/4]
    # x0_list = [wavelength/16, wavelength/8, wavelength*3/16]
    # alpha_list = [0.2, 0.4, 0.6, 0.8]
    alpha_list = [0.2, 0.4]
    t_list = np.linspace(0, 2*np.pi/Omega0*5, 50)
    fit = True

    # signal_1D = eval_signal_1D(t_list, alpha_list, x0_list)
    signal_3D = eval_signal_3D(t_list, alpha_list, x0_list)
    signal = signal_3D

    fig, axs = plt.subplots(len(alpha_list), len(x0_list), sharex=True, sharey=True, figsize=(3*len(x0_list), 2*len(alpha_list)), layout='tight')
    for i, alpha in enumerate(alpha_list):
        for j, x0 in enumerate(x0_list):
            axs[i, j].plot(t_list, signal[i, j], 'o', markersize=3)
            if fit:
                popt, pcov = curve_fit(damped_oscillation, t_list, signal[i, j], p0=[Omega0, 2*np.pi/Omega0*3])
                t_list_2 = np.linspace(t_list[0], t_list[-1], int(t_list[-1]/(2*np.pi/Omega0)*100))
                axs[i, j].plot(t_list_2, damped_oscillation(t_list_2, *popt))
                tau = popt[1]
                if tau < t_list[-1] * 3:
                    # otherwise the fit won't be accurate, since thre's not enough amplitude change within the fitted time range
                    axs[i, j].text(0.75, 0.9, r'$\Omega_0\tau=2\pi\times${:.2f}'.format(Omega0*popt[1]/(2*np.pi)), horizontalalignment='center', verticalalignment='center', transform=axs[i, j].transAxes)

            if i == 0:
                axs[i, j].set_title('x0 = {:.3f} $\lambda$'.format(x0/wavelength))

            if i == len(alpha_list) - 1:
                axs[i, j].set_xlabel('Time (ms)')

            if j == 0:
                axs[i, j].set_ylabel('Integrated signal\n'+r'($\alpha$ = {:.1f})'.format(alpha))            

    plt.show()