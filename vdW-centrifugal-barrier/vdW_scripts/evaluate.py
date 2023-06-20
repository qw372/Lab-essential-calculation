import numpy as np
from sympy.physics.wigner import clebsch_gordan, wigner_3j, wigner_6j, wigner_9j
from scipy import integrate

from .constants import *

__all__ = ['sort_eigenstates', 'unitarity_rate']

    
def sort_eigenstates(energies: np.ndarray, states: np.ndarray, assert_unique: bool = False) -> tuple[np.ndarray, np.ndarray]:
    ''' 
    Sort states to remove false avoided crossings.
    Adapted from https://github.com/qw372/Diatomic-Py-SrF/blob/master/diatomic/evaluation.py

    This is a function to ensure that all eigenstates plotted change
    adiabatically, it does this by assuming that step to step the eigenstates
    should vary by only a small amount (i.e. that the  step size is fine) and
    arranging states to maximise the overlap one step to the next.

    Args:
        energies (np.ndarray): 2D np array containing the eigenenergies, each row is a 1D array of eigenenergies under one condition
        states (np.ndarray): 3D np array containing the states, states[x, :, i] corresponds to energies[x, i]
        assert_unique (bool): If True, assert that there are no duplicate eigenstates (default False)
    Returns:
        energies (np.ndarray): 2D np array containing the eigenenergies, each row is a 1D array of eigenenergies under one condition but the order in each row is sorted
        states (np.ndarray): 3D np array containing the states, in the same order as energies E[x, j] -> States[x, :, j]
    '''

    assert len(energies.shape) == 2 # assert energies is 1D
    assert len(states.shape) == 3 # assert states is 2D
    assert energies.shape[0] == states.shape[0] # assert the number of energies is the same as number of states
    assert energies.shape[1] == states.shape[2] # assert the number of energies is the same as number of states

    number_iterations = energies.shape[0]

    # This loop sorts the eigenstates such that they maintain some continuity. 
    # Each eigenstate should be chosen to maximise the overlap with the previous.
    for i in range(number_iterations-1, 0, -1):

        #calculate the overlap of the ith and jth eigenstates
        overlaps = np.einsum('ij,ik->jk', np.conjugate(states[i, :, :]), states[i-1, :, :])
        orig_states = states[i-1, :, :].copy()
        orig_energies = energies[i-1, :].copy()

        ls = np.argmax(np.abs(overlaps), axis=1) # save the location of maximums into array ls

        if assert_unique:
            assert np.unique(ls).shape == ls.shape # assert there's no duplicates in ls

        for k, l in enumerate(ls):
            if l!=k:
                energies[i-1, k] = orig_energies[l].copy()
                states[i-1, :, k] = orig_states[:, l].copy()

    return (energies, states)

def unitarity_rate(T: float, centrifugal_barrier: float) -> float:
    """
    Calculate the unitarity collision rate a molecular sample at temperature T.
    See Thomas's writeup.

    Args:
        T: molecule temperature in uK
        centrifugal_barrier: centrifugal barrier in uK
    Returns:
        unitarity rate in cm^3/s
    """
    
    assert T > 0 and centrifugal_barrier >= 0

    centrifugal_barrier /= 1e6 # convert from uK to K
    T /= 1e6 # convert from uK to K

    hbar = 1.054571817e-34 # SI units, reduced Planck constant
    mu_SI = mu*9.1093837e-31 # SI units, reduced mass of SrF
    kB = 1.380649e-23 # SI units, Boltzmann constant

    factor = np.sqrt(2/np.pi) * (4*np.pi*hbar**2/mu_SI**2) * (mu_SI/kB/T)**(3/2)
    vc = np.sqrt(2*kB*centrifugal_barrier/mu_SI)

    integrand = lambda vr: vr * np.exp(-mu_SI*vr**2/2/kB/T)
    integral = integrate.quad(integrand, vc, np.inf)[0]

    rate = factor * integral * 1e6 # convert from m^3/s to cm^3/s

    return rate


# print(unitarity_rate(40, 0))