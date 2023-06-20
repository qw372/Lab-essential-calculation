import numpy as np
from sympy.physics.wigner import clebsch_gordan, wigner_3j, wigner_6j, wigner_9j
from scipy import integrate
from .constants import *


__all__ = ['centrifugal_coefficient_matrix', 'centrifugal_energy_matrix', 'centrigufal_barrier']


def centrifugal_coefficient_matrix(JA: int, JB: int, Lmax: int) -> float:
    """
    Args:
        JA: molecule A's rotational angular momentum
        JB: molecule B's rotational angular momentum
        Lmax: max partial wave to include
    Returns:
        centrifugal energy coefficient in atomic units
    """

    length = int((2*JA+1) * (2*JB+1) * (Lmax+1)**2)
    C = np.zeros((length, length), dtype=np.float64)
    i = 0
    for MA in np.arange(-JA, JA+1):
        for MB in np.arange(-JB, JB+1):
            for L in np.arange(0, Lmax+1):
                for ML in np.arange(-L, L+1):
                    C[i, i] = L*(L+1) / (2*mu)

                    i += 1

    return C

def centrifugal_energy_matrix(JA: int, JB: int, Lmax: int, R: float, unit_in_uK: bool = False) -> float:
    """
    Args:
        JA: molecule A's rotational angular momentum
        JB: molecule B's rotational angular momentum
        Lmax: max partial wave to include
        R: intermolecular distance in atomic units
        unit_in_uK: convert output's unit to uK
    Returns:
        centrifugal energy in atomic units (or uK)
    """

    C = centrifugal_coefficient_matrix(JA, JB, Lmax) / R**2

    return C * Hartree_to_uk if unit_in_uK else C


def centrigufal_barrier(C6: float, l: float, unit_in_uK: bool = False) -> float:
    """
    Calculate centrifugal barrier under a well-defined van der Waals interaction V(R) = -C6/R^6 + l(l+1)/(2*mu*R^2)

    Args:
        C6: van der Waals coefficient in atomic units
        l: partial wave angular momentum
        unit_in_uK: convert output's unit to uK
    Returns:
        centrifugal barrier in atomic units
    """    

    assert l >= 0 and C6 >= 0

    b = (l*(l+1)/mu)**1.5 / np.sqrt(54*C6)
 
    return b * Hartree_to_uk if unit_in_uK else b

    
# C = centrifugal_coefficient_matrix(1, 1, 5)
# print(C)
