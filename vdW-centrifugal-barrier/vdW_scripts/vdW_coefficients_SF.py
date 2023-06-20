import numpy as np
from sympy.physics.wigner import clebsch_gordan, wigner_3j, wigner_6j, wigner_9j
from scipy import integrate

from .constants import *

__all__ = ['W_SF_coefficient', 'W_SF', 'W_SF_coefficient_matrix', 'W_SF_matrix']


"""
Calculate second-order perturbation matrix element for van der Waals interaction between two N=1 SrF molecules (in space-fixed frame).

Following Eq. (82) from https://arxiv.org/pdf/1703.02833.pdf
"""


def _reduced_EDM(J: int, Omega: int, Jprime: int, Omegaprime: int, d: float) -> float:
    """
    Calculate reduced EDM matrix element <J, Omage1 || d || Jprime, Omegaprime> in atomic units

    Args:
        J, Jprime: rotational angular momentum
        Omega, Omegaprime: projection of rotational angular momentum on internuclear axis
        d: body-frame dipole moment in atomic units
    Returns:
        reduced EDM matrix element in atomic units
    """

    return (-1)**(J-Omega) * np.sqrt((2*J+1)*(2*Jprime+1)) * wigner_3j(J, 1, Jprime, -Omega, 0, Omegaprime) * d

def W_SF_coefficient(JA: int, MA: int, MA_prime: int, JB: int, MB: int, MB_prime: int, 
                     L: int, ML: int, L_prime: int, ML_prime: int) -> float:
    """
    Calculate second-order perturbation theory matrix element for van der Waals interaction coefficient C6
    between two N=1 SrF molecules, i.e., <JA, MA_prime, JB, MB_prime, L, L_prime | W | JA, MA, JB, MB, L, L>.

    Following Eq. (82) from https://arxiv.org/pdf/1703.02833.pdf

    Args:
        JA, MA, MA_prime: molecule A's rotational angular momentum
        JB, MB, MB_prime: molecule B's rotational angular momentum
        L, ML, L_prime, ML_prime: partial wave angular momentum
    Returns:
        W: matrix element of van der Waals interaction coefficient C6 in atomic units
    """

    assert JA >= 0 and JB >= 0 and L >= 0 and L_prime >= 0
    assert np.abs(MA) <= JA and np.abs(MB) <= JB and np.abs(ML) <= L and np.abs(ML_prime) <= L_prime

    if MA_prime + MB_prime + ML_prime != MA + MB + ML:
        return 0

    W = -1
    W *= (-1)**(lB_prime + lB + 2*JA + 2*JB)
    W *= np.sqrt(np.math.factorial(2*l_prime+1) * np.math.factorial(2*l+1)) 
    W /= np.sqrt(np.math.factorial(2*lA_prime) * np.math.factorial(2*lB_prime))
    W /= np.sqrt(np.math.factorial(2*lA) * np.math.factorial(2*lB))
    
    W_partial = 0
    qA = MA_prime - MA
    qB = MB_prime - MB
    q = qA + qB
    for kA in np.arange(np.abs(lA-lA_prime), lA+lA_prime+1):
        if np.abs(qA) > kA:
            continue

        for kB in np.arange(np.abs(lB-lB_prime), lB+lB_prime+1):
            if np.abs(qB) > kB:
                continue

            for k in np.arange(np.abs(kA-kB), kA+kB+1):
                if np.abs(q) > k:
                    continue

                factor = clebsch_gordan(l_prime, l, k, 0, 0, 0) 
                factor *= np.sqrt((2*L_prime+1)/(2*L+1)) * clebsch_gordan(L_prime, k, L, 0, 0, 0) * clebsch_gordan(L_prime, k, L, ML_prime, q, ML)
                factor *= ((-1)**(kA+kB)) * (2*kA+1) * (2*kB+1) * clebsch_gordan(kA, kB, k, qA, qB, q)
                if factor == 0:
                    continue

                s = 0
                for JA_doubleprime in np.arange(np.abs(JA-lA), JA+lA+1):
                    if JA == JA_doubleprime:
                        continue

                    for JB_doubleprime in np.arange(np.abs(JB-lB), JB+lB+1):
                        if JB == JB_doubleprime:
                            continue

                        Delta = JA_doubleprime*(JA_doubleprime+1)*B_rot - JA*(JA+1)*B_rot + JB_doubleprime*(JB_doubleprime+1)*B_rot - JB*(JB+1)*B_rot
                        assert Delta != 0

                        s_partial = 1

                        s_partial *= _reduced_EDM(JA, 0, JA_doubleprime, 0, d_EDM) * _reduced_EDM(JA_doubleprime, 0, JA, 0, d_EDM)
                        s_partial *= _reduced_EDM(JB, 0, JB_doubleprime, 0, d_EDM) * _reduced_EDM(JB_doubleprime, 0, JB, 0, d_EDM)
                        s_partial /= Delta
                        s_partial *= wigner_9j(lA_prime, lA, kA, lB_prime, lB, kB, l_prime, l, k)
                        s_partial *= wigner_6j(lA, lA_prime, kA, JA, JA, JA_doubleprime)
                        s_partial *= wigner_6j(lB, lB_prime, kB, JB, JB, JB_doubleprime)
                        s_partial *= clebsch_gordan(JA, kA, JA, MA, qA, MA_prime) * clebsch_gordan(JB, kB, JB, MB, qB, MB_prime)
                        s_partial /= np.sqrt((2*JA+1)*(2*JB+1))

                        s += s_partial

                W_partial += factor * s
                # print(s)

    W *= W_partial

    return W

def W_SF(JA: int, MA: int, MA_prime: int, JB: int, MB: int, MB_prime: int, 
         L: int, ML: int, L_prime: int, ML_prime: int, R: float) -> float:
    """
    Calculate second-order perturbation theory matrix element for van der Waals interaction
    between two N=1 SrF molecules, i.e., <JA, MA_prime, JB, MB_prime, L, L_prime | W | JA, MA, JB, MB, L, L>.

    Following Eq. (82) from https://arxiv.org/pdf/1703.02833.pdf

    Args:
        JA, MA, MA_prime: molecule A's rotational angular momentum
        JB, MB, MB_prime: molecule B's rotational angular momentum
        L, ML, L_prime, ML_prime: partial wave angular momentum
        R: intermolecular distance in atomic units
    Returns:
        W: matrix element of van der Waals interaction in atomic units
    """

    return W_SF_coefficient(JA, MA, MA_prime, JB, MB, MB_prime, L, ML, L_prime, ML_prime) / R**6

def W_SF_coefficient_matrix(JA: int, JB: int, Lmax: int) -> np.ndarray:
    """
    Calculate second-order perturbation theory matrix for van der Waals interaction coefficient C6
    between two N=1 SrF molecules, i.e., <JA, MA_prime, JB, MB_prime, L, L_prime | W | JA, MA, JB, MB, L, L>.

    Args:
        JA: molecule A's rotational angular momentum
        JB: molecule B's rotational angular momentum
        Lmax: max partial wave to include
    Returns:
        W: matrix of van der Waals interaction coefficient C6 in atomic units
    """

    assert JA >= 0 and JB >= 0 and Lmax >= 0
    
    length = int((2*JA+1) * (2*JB+1) * (Lmax+1)**2)
    index_list = np.empty((length, 4), dtype=np.float64)
    i = 0
    for MA in np.arange(-JA, JA+1):
        for MB in np.arange(-JB, JB+1):
            for L in np.arange(0, Lmax+1):
                for ML in np.arange(-L, L+1):
                    index_list[i, 0] = MA
                    index_list[i, 1] = MB
                    index_list[i, 2] = L
                    index_list[i, 3] = ML

                    i += 1

    W = np.zeros((length, length), dtype=np.float64)
    for i, index in enumerate(index_list):
        MA = index[0]
        MB = index[1]
        L = index[2]
        ML = index[3]

        W[i, i] = W_SF_coefficient(JA, MA, MA, JB, MB, MB, L, ML, L, ML) / 2 # divide by 2 here becasue of W += W.T step later

        for j in np.arange(i+1, length):
            MA_prime = index_list[j, 0]
            MB_prime = index_list[j, 1]
            L_prime = index_list[j, 2]
            ML_prime = index_list[j, 3]

            W[i, j] = W_SF_coefficient(JA, MA, MA_prime, JB, MB, MB_prime, L, ML, L_prime, ML_prime)

        # print(i)

    W += W.T

    return W

def W_SF_matrix(JA: int, JB: int, Lmax: int, R: float, unit_in_uK: bool = False) -> np.ndarray:
    """
    Calculate second-order perturbation theory matrix for van der Waals interaction
    between two N=1 SrF molecules, i.e., <JA, MA_prime, JB, MB_prime, L, L_prime | W | JA, MA, JB, MB, L, L>.

    Args:
        JA: molecule A's rotational angular momentum
        JB: molecule B's rotational angular momentum
        Lmax: max partial wave to include
        R: intermolecular distance in atomic units
        unit_in_uK: convert output's unit to uK
    Returns:
        W: matrix of van der Waals interaction in atomic units (or uK)
    """

    W = W_SF_coefficient_matrix(JA, JB, Lmax) / R**6

    return W * Hartree_to_uk if unit_in_uK else W


# W = W_SF_coefficient_matrix(1, 1, 2)
# print(W)
