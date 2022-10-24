import numpy as np
from sympy.physics.wigner import *

Flist = [1, 2]
mlist = [-2, -1, 0, 1, 2]
qlist = [-1, 0, 1]

Fprime = 3 # excited state F
mprime = 3 # excited state m
Jprime = 3/2 # excited state J
J = 1/2 # ground state J
I = 3/2 # nuclear spin

s = 0
for F in Flist:
    for m in mlist:
        for q in qlist:
            s += (np.abs(wigner_3j(Fprime, 1, F, -mprime, q, m))**2*(2*Fprime+1)*(2*F+1)
                    *(np.abs(wigner_6j(Jprime, Fprime, I, F, J, 1)))**2)

print(s)
gamma = 2*np.pi*6.065e6 # in Hz
alpha = 1/137
c = 299792458
omega = 2*np.pi*c/780e-9 # in Hz, transition frequency

red_mat_elem_J = (gamma/(4/3*alpha*omega**3/c**2*s))**0.5
# print(red_mat_elem_J)

red_mat_elem_F = float(red_mat_elem_J*(2*F+1)**0.5*(2*Fprime+1)**0.5*wigner_6j(Jprime, Fprime, I, F, J, 1))
print(red_mat_elem_F)

I = 10e-6/(np.pi*0.01**2) # in W/m^2, laser intensity
epsilon = 8.8541878128e-12
E = (2*I/c/epsilon)**0.5
print(E)
hbar = 1.0545718e-34
electron = 1.60217662e-19 # in C, electron charge
rabifreq = E*red_mat_elem_F*electron/(hbar*gamma)
print(rabifreq)
