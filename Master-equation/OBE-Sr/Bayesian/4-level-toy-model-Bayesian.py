import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from multiprocessing import Pool
from skopt import gp_minimize
from skopt.plots import plot_convergence
import time


def masterequation(t, den, Ham, Gamma_21, Gamma_32, Gamma_34):
    """
    Master eqaution

    t: time
    den: density matrix (in rotating frame)
    Ham: Hamiltonian (in rotating frame)
    Gamma_21: decay rate from state 2 to state 1 
    Gamma_32: decay rate from state 3 to state 2 
    Gamma_34: decay rate from state 3 to state 4 
    """

    den = den.reshape((4,4))


    re = -1j*(np.matmul(Ham, den)-np.matmul(den, Ham))

    # Lindblad_21
    re += Gamma_21*np.array([[den[1, 1],  -den[0, 1]/2,  0,             0],
                             [-den[1, 0]/2, -den[1, 1],    -den[1, 2]/2,  -den[1, 3]/2],
                             [0,            -den[2, 1]/2,  0,             0],
                             [0,            -den[3, 1]/2,  0,             0]])
                    
    # Lindblad_32
    re += Gamma_32*np.array([[0,           0,             -den[0, 2]/2,  0],
                             [0,             den[2, 2],     -den[1, 2]/2,  0],
                             [-den[2, 0]/2,  -den[2, 1]/2,  -den[2, 2],    -den[2, 3]/2],
                             [0,             0,             -den[3, 2]/2,  0]])

    # Lindblad_34
    re += Gamma_34*np.array([[0,           0,             -den[0, 2]/2,  0],
                             [0,             0,             -den[1, 2]/2,  0],
                             [-den[2, 0]/2,  -den[2, 1]/2,  -den[2, 2],    -den[2, 3]/2],
                             [0,             0,             -den[3, 2]/2,  den[2, 2]]])

    return re.reshape(16)   

def solve_me(Omega_12, Delta_12):
    """
    Solve master equation
    Omega_12: Rabi freq for 1-to-2 transition
    Delta_12: detuning for 1-to-2 transition
    """

    Gamma_21 = 2*np.pi*32 # MHz, decay rate
    Gamma_32 = 2*np.pi*3 # MHz, decay rate
    Gamma_34 = 2*np.pi*0.5 # MHz, decay rate

    # artificially make them dependent to test Bayesian optimization
    Omega_23 = 20 - Omega_12 # MHz, Rabi freq
    Delta_23 = -Delta_12 # MHz, detuning

    if Omega_12 < 0:
        print(f"Omega_12 ({Omega_12}) less than 0.")
        return -1
    elif Omega_23 < 0:
        print(f"Omega_23 ({Omega_23}) less than 0.")
        return -1

    Ham = np.array([[0,          Omega_12/2,  0,                   0],
                    [Omega_12/2,  -Delta_12,   Omega_23/2,          0],
                    [0,           Omega_23/2,  -Delta_12-Delta_23,  0],
                    [0,           0,           0,                   0]])
    
    den_0 = np.array([[1+0j, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0]]) # initial density matrix
    den_0 = den_0.reshape(16)

    t_span = [0, 1000]

    sol = solve_ivp(lambda t, den: masterequation(t, den, Ham, Gamma_21, Gamma_32, Gamma_34), t_span=t_span, y0=den_0, vectorized=True)
    t = sol.t
    # pop_1 = np.abs(sol.y[0])**2 # population on state 1
    # pop_2 = np.abs(sol.y[5])**2 # population on state 2
    # pop_3 = np.abs(sol.y[10])**2 # population on state 3
    pop_4 = np.abs(sol.y[15])**2 # population on state 4

    return pop_4[-1]

def solve_me_for_bayesian(param):
    """
    Re-write solve_me function for bayesian optimization
    param: [Omega_12, Delta_12]
    """
    
    return -solve_me(param[0], param[1])

def solve_me_multiprocess(Omega_12_list, Delta_12_list):
    """
    Solve master euqation parallelly for various paraemters
    Omega_12_list: list of Omega_12 to solve
    Delta_12_list: list of Delta_12 to sovle
    """
    
    Omega_12_list2 = np.repeat(Omega_12_list, len(Delta_12_list))
    Delta_12_list2 = np.tile(Delta_12_list, len(Omega_12_list))
    param_list = list(zip(Omega_12_list2, Delta_12_list2))

    if __name__ == "__main__":
        with Pool(8) as p:
            result = np.array(p.starmap(solve_me, param_list))
            # print(result)

            pop_4_list = np.reshape(result, (-1, len(Delta_12_list)))

    return pop_4_list

def solve_me_bayesian(Omega_12_range, Delta_12_range):
    """
    Optimize parameters of master equation using Bayesian optimization
    Omega_12_range: range of Omega_12
    Delta_12_range: range of Delta_12

    Further optimization: parallelize this?
    """
    
    res = gp_minimize(solve_me_for_bayesian,                  # the function to minimize
                      [Omega_12_range, Delta_12_range],      # the bounds on each dimension of x
                      acq_func="EI",      # the acquisition function
                      n_calls=10,         # the number of evaluations of f
                      n_random_starts=1,  # the number of random initialization points
                      random_state=1234)   # the random seed

    return res

if __name__ == "__main__":
    Omega_12_range = (0.0, 13.0)
    Delta_12_range = (-10.0, 2.0)
    Omega_12_list = np.linspace(Omega_12_range[0], Omega_12_range[1], 20)
    Delta_12_list = np.linspace(Delta_12_range[0], Delta_12_range[1], 20)
    BayesianOptim = True

    if BayesianOptim:
        t = time.time()
        res = solve_me_bayesian(Omega_12_range, Delta_12_range)
        print("Running time: {:.2f} s.".format(time.time()-t))
        ax = plot_convergence(res)
        textstr = "Optimum param: ({:.2f}, {:.2f})\nOptimized func value: {:.2f}".format(res.x[0], res.x[1], res.fun)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.98, 0.97, textstr, transform=ax.transAxes, verticalalignment='top', horizontalalignment='right', bbox=props)
        plt.savefig("Bayesian-convergence-plot.png")
        plt.show()

    else:
        t = time.time()
        pop_4_list = solve_me_multiprocess(Omega_12_list, Delta_12_list)
        print("Running time: {:.2f} s.".format(time.time()-t))

        plt.imshow(pop_4_list[::-1], aspect='auto', extent=[Delta_12_list[0],Delta_12_list[-1],Omega_12_list[0],Omega_12_list[-1]])
        plt.xlabel(r"$\Delta_{12}$", fontsize=14)
        plt.ylabel(r"$\Omega_{12}$", fontsize=14)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('State 4 population', rotation=90, fontsize=12, labelpad=12)
        plt.savefig("2D-map.png")
        plt.show()