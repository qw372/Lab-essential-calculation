import numpy as np
from classes import *


LaserPower_W = 5 # laser power in W
LaserRadius_inve2_um2 = 40 # 1/e^2 radius of laser beam, in um
LaserIntensity_kW_invcm2 = (LaserPower_W/2/np.pi/(LaserRadius_inve2_um2/2)**2)/1e3*1e8 # peak Laser intensity in kW/cm^2
# print(LaserIntensity_kW_invcm2)

SrF_state_1 = Hunds_case_b_state(label="SrF XSigma", Lambda=0, N=1, S=1/2, J=3/2, I=1/2, F=2)
SrF_state_2 = Hunds_case_b_state(label="SrF XSigma", Lambda=0, N=1, S=1/2, J=3/2, I=1/2, F=1, JMixing=[{"mixing coefficient":0.8880, "J":3/2}, {"mixing coefficient":0.45984345162, "J":1/2}])
SrF_state_3 = Hunds_case_b_state(label="SrF XSigma", Lambda=0, N=1, S=1/2, J=1/2, I=1/2, F=0)
SrF_state_4 = Hunds_case_b_state(label="SrF XSigma", Lambda=0, N=1, S=1/2, J=1/2, I=1/2, F=1, JMixing=[{"mixing coefficient":-0.45984345162, "J":3/2}, {"mixing coefficient":0.8880, "J":1/2}])

SrF_state_5 = Hunds_case_b_state(label="SrF XSigma", Lambda=0, N=0, S=1/2, J=1/2, I=1/2, F=1)
SrF_state_6 = Hunds_case_b_state(label="SrF XSigma", Lambda=0, N=0, S=1/2, J=1/2, I=1/2, F=0)

for state in [SrF_state_1, SrF_state_2, SrF_state_3, SrF_state_4, SrF_state_5, SrF_state_6]:
    SrFacStarkShift(state, LaserWavelength_nm=1064, LaserIntensity_kW_invcm2=LaserIntensity_kW_invcm2, print_dipole_moment=False, print_polarizability=True, print_stark_shift=True)

Rb_state_1 = Rb_ground_state(label="Rb 5S", energy=0, J=1/2, I=3/2, F=2)
Rb_state_2 = Rb_ground_state(label="Rb 5S", energy=0, J=1/2, I=3/2, F=1)

for state in [Rb_state_1, Rb_state_2]:
    RbacStarkShift(state, LaserWavelength_nm=1064, LaserIntensity_kW_invcm2=LaserIntensity_kW_invcm2, print_polarizability=True, print_stark_shift=True)