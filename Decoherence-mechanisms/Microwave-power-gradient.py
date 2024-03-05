import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

wavelength = 20 # mm, microwave wavelength
k = 2*np.pi/wavelength
omega0 = 2*np.pi*50 # kHz, peak Rabi frequency
sigma = 0.15 # mm, ODT cloud size

def rabi_freq(x, x0, alpha):
    # return omega0 * np.sqrt(1 + alpha**2 + 2 * alpha * np.cos(2 * k * (x - x0))) / (1 + alpha)
    return omega0 * (1 + alpha * np.sin(2 * k * (x-x0)))

def Rabi_oscillation(t, x, x0, alpha):
    return 0.5 + 0.5 * np.cos(rabi_freq(x, x0, alpha) * t)

def gaussian_dist(x):
    return np.exp(-0.5 * (x/sigma)**2) / (sigma * np.sqrt(2*np.pi))

def integrand(t, x, x0, alpha):
    return Rabi_oscillation(t, x, x0, alpha) * gaussian_dist(x)

def result(t_list, x0, alpha):
    result = np.empty(len(t_list))
    for i, t in enumerate(t_list):
        result[i] = integrate.quad(lambda x: integrand(t=t, x=x, x0=x0, alpha=alpha), -np.Infinity, np.Infinity)[0]

    return result

x0_list = [0, wavelength/16, wavelength/8, wavelength/4]
alpha_list = [0.2, 0.4, 0.6, 0.8]
t_list = np.linspace(0, 2*np.pi/omega0*10, 200)

fig, axs = plt.subplots(len(alpha_list), len(x0_list), sharex=True, sharey=True, figsize=(12, 8), layout='tight')

for i, alpha in enumerate(alpha_list):
    for j, x0 in enumerate(x0_list):
        axs[i, j].plot(t_list, result(t_list, x0, alpha), label=r'$\alpha$ = {:.1f}'.format(alpha))

        if i == 0:
            axs[i, j].set_title('x0 = {:.3f} $\lambda$'.format(x0/wavelength))

        if i == len(alpha_list) - 1:
            axs[i, j].set_xlabel('Time (ms)')

        if j == 0:
            axs[i, j].set_ylabel('Integrated signal')

        if j == len(x0_list) - 1:
            axs[i, j].legend(loc='upper right')

plt.legend()
plt.show()