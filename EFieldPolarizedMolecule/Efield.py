import numpy as np
from numpy import linalg
from sympy.physics.wigner import *
import matplotlib.pyplot as plt

N_to_calculate = 100
N_to_plot = 5
E_field_list = np.linspace(0, 200, 200)
Ham_diag = np.diag([n*(n+1) for n in range(N_to_calculate)])
off_diag_array = np.array([float(wigner_3j(n, 1, n+1, 0, 0, 0)**2)*np.sqrt(2*n+1)*np.sqrt(2*n+3) for n in range(N_to_calculate-1)])
dipole_moment_mat = np.diag(off_diag_array, k=1) + np.diag(off_diag_array, k=-1)

eigenvalue_list = np.array([])
dipole_moment_list = np.array([])
for E in E_field_list:
    Ham = Ham_diag + dipole_moment_mat*E

    eigenvalues, eigenstates = linalg.eigh(Ham)
    eigenvalues = eigenvalues[:N_to_plot]
    eigenstates = eigenstates.T[:N_to_plot, ::]
    eigenvalue_list = np.append(eigenvalue_list, eigenvalues)
    for state in eigenstates:
        dipole_moment_list = np.append(dipole_moment_list, np.matmul(np.matmul(state.T, dipole_moment_mat), state))

eigenvalue_list = eigenvalue_list.reshape((len(E_field_list), -1)).T
dipole_moment_list = dipole_moment_list.reshape(len(E_field_list), -1).T

ax1 = plt.subplot(1, 2, 1)
for eigenvalue_sublist in eigenvalue_list:
    ax1.plot(E_field_list, eigenvalue_sublist)
ax1.set_xlabel("Electric field / (B/d)")
ax1.set_ylabel("Energy / B")
ax1.grid()

ax2 = plt.subplot(1, 2, 2)
for dipole_moment_sublist in dipole_moment_list:
    ax2.plot(E_field_list, dipole_moment_sublist)
ax2.set_xlabel("Electric field / (B/d)")
ax2.set_ylabel("Dipole moment / d")
ax2.grid()

plt.show()