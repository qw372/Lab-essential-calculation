__all__ = ['lA', 'lB', 'l', 'lA_prime', 'lB_prime', 'l_prime', 'B_rot', 'd_EDM', 'mu', 'Hartree_to_uk']

lA = 1 # index for (electrostatic) multipole expansion of molecule A
lB = 1 # index for (electrostatic) multipole expansion of molecule B
l = lA + lB
lA_prime = 1 # index for (electrostatic) multipole expansion of molecule A
lB_prime = 1 # index for (electrostatic) multipole expansion of molecule B
l_prime = lA_prime + lB_prime

B_rot = 7.60274/6579.683920502e3 # rotational constant of SrF in atomic units
d_EDM = 3.4963/2.541746473 # body-frame dipole moment of SrF in atomic units
mu = 53.45*1822.89 # reduced mass of SrF in atomic units
Hartree_to_uk = 315775.02480407e6 # factor to convert Hartree to uK