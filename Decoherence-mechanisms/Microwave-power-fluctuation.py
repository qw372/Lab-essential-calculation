import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

Omega0 = 2*np.pi*10e3 # Hz, center Rabi frequency
sigma_Omega = Omega0 * 0.3 # std of Rabi frequency

cycle = 5
t_span = [0, 2*np.pi/Omega0*cycle] # s, time span
t_eval = np.linspace(t_span[0], t_span[1], cycle*20)
# t_noise = np.linspace(t_span[0], t_span[1], cycle*100)

rng = np.random.default_rng(seed=12345)

def evaluete():
    # Omega_noise = rng.normal(Omega0, sigma_Omega, len(t_noise))

    def func(t, y):
        # Omega = np.interp(t, t_noise, Omega_noise)
        Omega = rng.normal(Omega0, sigma_Omega)
        return [1/2j*Omega*y[1], 1/2j*Omega*y[0]]

    sol = solve_ivp(func, t_span=t_span, y0=[0+0j, 1+0j], t_eval=t_eval)
    return np.abs(sol.y[1])**2

p = np.zeros(len(t_eval))
num = 10
for i in range(num):
    if i == 0:
        p = evaluete()
    else:
        p = np.vstack((p, evaluete()))

plt.errorbar(t_eval, np.mean(p, axis=0), yerr=np.std(p, axis=0)/np.sqrt(num))
# plt.plot(Omega_noise+Omega0)
plt.show()