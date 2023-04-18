import numpy as np
from sympy import sqrt
from sympy.physics.wigner import *

d = 1.3755
LaserWavelength_nm = 1064
nmToAtomicUnit = 1/0.0529177210903
SpeedOfLight = 137.035999084
LaserWavelength = LaserWavelength_nm*nmToAtomicUnit
hbar = 1
LaserEnergy = SpeedOfLight/LaserWavelength*hbar*2*np.pi

Lambda = 0
N = 1
S = 1/2
J = 1/2
I = 1/2
F = 1

Lambda_prime = 0
K = 1
partial_polarizability = 0
for N_prime in np.arange(np.abs(Lambda_prime), N+1+1):
    for J_prime in np.arange(np.abs(N_prime-S), N_prime+S+1):
        for F_prime in np.arange(np.abs(J_prime-I), J_prime+I+1):
            if (np.abs(N_prime-N)<0.1) and (np.abs(J_prime-J)<0.1) and (np.abs(F_prime-F)<0.1):
                continue

            partial_dipole = 1
            partial_dipole *= ((-1)**(F+J_prime+1+I))*sqrt(2*F+1)*sqrt(2*F_prime+1)*wigner_6j(J, F, I, F_prime, J_prime, 1)
            partial_dipole *= ((-1)**(J+N_prime+1+S))*sqrt(2*J_prime+1)*sqrt(2*J+1)*wigner_6j(N, J, S, J_prime, N_prime, 1)
            partial_dipole *= ((-1)**(N_prime-Lambda_prime))*wigner_3j(N_prime, 1, N, -Lambda_prime, 0, Lambda)*sqrt(2*N_prime+1)*sqrt(2*N+1)

            print(partial_dipole**2)

            partial_polarizability += ((-1)**(F_prime+F))*wigner_6j(1, 1, K, F, F, F_prime)*(partial_dipole**2)

# print(partial_polarizability)
partial_polarizability *= sqrt(2*K+1)*((-1)**K)*(-2/LaserEnergy)
partial_polarizability *= sqrt(2*F/(F+1)/(2*F+1))
# partial_polarizability = float(partial_polarizability)

print("")
print(partial_polarizability)