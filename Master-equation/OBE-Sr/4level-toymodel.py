import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time

def Lindblad_21(den):

    Gamma_21 = 2*np.pi*32 # MHz, decay rate
    return Gamma_21*np.array([[den[1, 1],  -den[0, 1]/2,  0,             0],
                            [-den[1, 0]/2, -den[1, 1],    -den[1, 2]/2,  -den[1, 3]/2],
                            [0,            -den[2, 1]/2,  0,             0],
                            [0,            -den[3, 1]/2,  0,             0]])

def Lindblad_32(den):

    Gamma_32 = 2*np.pi*3 # MHz, decay rate
    return Gamma_32*np.array([[0,           0,             -den[0, 2]/2,  0],
                            [0,             den[2, 2],     -den[1, 2]/2,  0],
                            [-den[2, 0]/2,  -den[2, 1]/2,  -den[2, 2],    -den[2, 3]/2],
                            [0,             0,             -den[3, 2]/2,  0]])

def Lindblad_34(den):

    Gamma_34 = 2*np.pi*0.5 # MHz, decay rate
    return Gamma_34*np.array([[0,           0,             -den[0, 2]/2,  0],
                            [0,             0,             -den[1, 2]/2,  0],
                            [-den[2, 0]/2,  -den[2, 1]/2,  -den[2, 2],    -den[2, 3]/2],
                            [0,             0,             -den[3, 2]/2,  den[2, 2]]])

def Hamiltonian():

    Omega_12 = 10 # MHz, Rabi freq
    Omega_23 = 3 # MHz, Rabi freq

    Delta_12 = 0 # MHz, detuning
    Delta_23 = 0 # MHz, detuning

    return np.array([[0,          Omega_12/2,  0,                   0],
                    [Omega_12/2,  -Delta_12,   Omega_23/2,          0],
                    [0,           Omega_23/2,  -Delta_12-Delta_23,  0],
                    [0,           0,           0,                   0]])


def masterequation(t, den):
    """
    t: time
    den: density matrix (in rotating frame)
    """

    den = den.reshape((4,4))

    re = -1j*(np.matmul(Hamiltonian(), den)-np.matmul(den, Hamiltonian()))
    re += Lindblad_21(den)
    re += Lindblad_32(den)
    re += Lindblad_34(den)

    return re.reshape(16)   


den_0 = np.array([[1+0j, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0]]) # initial density matrix
den_0 = den_0.reshape(16)
t_span = [0, 1000]

t0 = time.time()
sol = solve_ivp(masterequation, t_span=t_span, y0=den_0, t_eval=np.linspace(t_span[0], t_span[1], 200), vectorized=True)
print("evaluation time: {:.2f} s.".format(time.time()-t0))

t = sol.t
pop_1 = np.abs(sol.y[0])**2 # population on state 1
pop_2 = np.abs(sol.y[5])**2 # population on state 2
pop_3 = np.abs(sol.y[10])**2 # population on state 3
pop_4 = np.abs(sol.y[15])**2 # population on state 4

plt.plot(t, pop_1, label="state 1")
plt.plot(t, pop_2, label="state 2")
plt.plot(t, pop_3, label="state 3")
plt.plot(t, pop_4, label="state 4")
plt.legend()
plt.ylabel("population")
plt.xlabel("evaluation time / us")
plt.grid()
# plt.savefig("latest.png")
plt.show()