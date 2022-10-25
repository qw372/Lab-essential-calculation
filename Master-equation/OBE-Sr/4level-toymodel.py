import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def Lindblad_21(den):

    Gamma_21 = 1 # MHz, decay rate
    return Gamma_21*np.array([[den[1, 1],  -den[0, 1]/2,  0,             0],
                            [-den[1, 0]/2, -den[1, 1],    -den[1, 2]/2,  -den[1, 3]/2],
                            [0,            -den[2, 1]/2,  0,             0],
                            [0,            -den[3, 1]/2,  0,             0]])

def Lindblad_32(den):

    Gamma_32 = 1 # MHz, decay rate
    return Gamma_32*np.array([[0,           0,             -den[0, 2]/2,  0],
                            [0,             den[2, 2],     -den[1, 2]/2,  0],
                            [-den[2, 0]/2,  -den[2, 1]/2,  -den[2, 2],    -den[2, 3]/2],
                            [0,             0,             -den[3, 2]/2,  0]])

def Lindblad_34(den):

    Gamma_34 = 1 # MHz, decay rate
    return Gamma_34*np.array([[0,           0,             -den[0, 2]/2,  0],
                            [0,             0,             -den[1, 2]/2,  0],
                            [-den[2, 0]/2,  -den[2, 1]/2,  -den[2, 2],    -den[2, 3]/2],
                            [0,             0,             -den[3, 2]/2,  den[2, 2]]])

def Hamiltonian():

    Omega_12 = 1 # MHz, Rabi freq
    Omega_23 = 1 # MHz, Rabi freq

    Delta_12 = 1 # MHz, detuning
    Delta_23 = 1 # MHz, detuning

    omega_4 = 3e8 # energy of state |4>, mneasured from state |1>

    return np.array([[0,          Omega_12/2,  0,                   0],
                    [Omega_12/2,  -Delta_12,   Omega_23/2,          0],
                    [0,           Omega_23/2,  -Delta_12-Delta_23,  0],
                    [0,           0,           0,                   omega_4]])


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
t_span = [0, 10]

sol = solve_ivp(masterequation, t_span=t_span, y0=den_0, t_eval=np.linspace(t_span[0], t_span[1], 100), vectorized=True)
print(sol.t)
print(sol.y)