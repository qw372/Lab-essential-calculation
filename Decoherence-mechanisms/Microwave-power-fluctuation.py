import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

t_span = [0, 0.5e-3]
t_eval = np.linspace(t_span[0], t_span[1], 1000)
t_eval2 = np.linspace(t_span[0], t_span[1], 3)

rng = np.random.default_rng(seed=12345)
Omega0 = 2*np.pi*10e3
Omega1 = 2*np.pi*2e3
# Omega_noise = Omega1 * (rng.random(len(t_eval2)) - 0.5)
Omega_noise = Omega1 * np.array([1, 0.5, 0])
# Omega = lambda t: np.interp(t, t_eval2, Omega_noise) + Omega0

def func(t, y):
    Omega = np.interp(t, t_eval2, Omega_noise) + Omega0
    return [1/2j*Omega*y[1], 1/2j*Omega*y[0]]

sol = solve_ivp(func, t_span=t_span, y0=[0+0j, 1+0j], t_eval=t_eval)
# print(np.abs(sol.y[1])**2)

plt.plot(np.abs(sol.y[1])**2)
# plt.plot(Omega_noise+Omega0)
plt.show()