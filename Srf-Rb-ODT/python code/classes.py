import numpy as np
from sympy.physics.wigner import *

def wigner_eckart_coefficient(j, m, k, p, jprime, mprime):
    r = (-1)**(j-m)*wigner_3j(j, k, jprime, -m, p, mprime)
    return float(r)

def spectator_theorem_coefcient(j1, j2, j12, k1, j1prime, j2prime, j12prime):
    if np.abs(j2 - j2prime) < 1e-6:
        r = (-1)**(j12prime+j1+k1+j2)*np.sqrt(2*j12+1)*np.sqrt(2*j12prime+1)*wigner_6j(j1prime, j12prime, j2, j12, j1, k1)
        return float(r)
    else:
        return 0

class Hunds_case_a_state:
    def __init__(self, label, Lambda, S, Sigma, J, Omega, I, F, mF):
        self.label = label
        self.Lambda = Lambda
        self.S = np.abs(S)
        self.Sigma = Sigma
        self.J = np.abs(J)
        self.Omega = Omega
        self.I = np.abs(I)
        self.F = np.abs(F)
        self.mF = mF

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        pass

    def calculate_dipole_moment(self, dipole_moment_Lambda, state_prime):
        Lambda = self.Lambda
        S = self.S
        Sigma = self.Sigma
        J = self.J
        Omega = self.Omega
        I = self.I
        F = self.F

        Lambda_prime = state_prime.Lambda
        S_prime = state_prime.S
        Sigma_prime = state_prime.Sigma
        J_prime = state_prime.J
        Omega_prime = state_prime.Omega
        I_prime = state_prime.I
        F_prime = state_prime.F

        if np.abs(Sigma - Sigma_prime) > 1e-1:
            # print("Angular momentum Sigma doesn't match.")
            return 0

        if np.abs(S - S_prime) > 1e-1:
            # print("Angular momentum S don't match.")
            return 0

        if np.abs(I -I_prime) > 1e-1:
            # print("Angular momentum I don't match.")
            return 0

        dipole = dipole_moment_Lambda
        dipole *= (-1)**(F+2*J_prime+1+I-Omega_prime)
        dipole *= np.sqrt(2*F_prime+1)*np.sqrt(2*F+1)*float(wigner_6j(J, F, I, F_prime, J_prime, 1))
        dipole *= np.sqrt(2*J_prime+1)*np.sqrt(2*J+1)*float(wigner_3j(J_prime, 1, J, -Omega_prime, -Omega+Omega_prime, Omega))

        return dipole


class Hunds_case_b_state:
    def __init__(self, label, Lambda, N, S, J, I, F, mF, JMixing=[]):
        self.label = label
        self.Lambda = Lambda
        self.N = np.abs(N)
        self.S = np.abs(S)
        self.J = np.abs(J)
        self.I = np.abs(I)
        self.F = np.abs(F)
        self.mF = mF

        self.convert_to_Hunds_case_a(JMixing)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        pass

    def convert_to_Hunds_case_a(self, JMixing):
        # JMixing is a list of dictionaries that saves J mixing infomation,
        # in form of [{"mixing coefficient": a number, "J": J}, ...]

        # define  list of dictionaries that saves Hund's a states mixd in this Hund's b state, 
        # in form of [{"coefficient": a number, "state": a Hund's a state}, ...]
        self.Hunds_case_a_list = []

        if JMixing:
            for mixed_b_state_info in JMixing:
                # mixed_b_state_info is in form of {"mixing coefficient": a number, "J": J}

                b_state = Hunds_case_b_state(label=self.label, Lambda=self.Lambda, N=self.N, S=self.S, J=mixed_b_state_info["J"], I=self.I, F=self.F, mF=self.mF)
                for a_state in b_state.Hunds_case_a_list:
                    # a_state is in form of {"coefficient": a number, "state": a Hund's a state}

                    d = {}
                    d["coefficient"] = a_state["coefficient"]*mixed_b_state_info["mixing coefficient"]
                    d["state"] = a_state["state"]
                    self.Hunds_case_a_list.append(d)

        else:
            # follow Eq. 2.1 of Eric's thesis 
            # https://cpb-us-w2.wpmucdn.com/voices.uchicago.edu/dist/2/2979/files/2020/12/norrgardthesis2016_04_26.pdf

            for Sigma in np.arange(-self.S, self.S+1):
                coefficient = np.sqrt(2*self.N+1)*float(wigner_3j(self.S, self.N, self.J, Sigma, self.Lambda, -Sigma-self.Lambda))*((-1)**(self.J+Sigma+self.Lambda))
                if coefficient != 0:
                    # print(coefficient)
                    d = {}
                    d["coefficient"] = coefficient
                    d["state"] = Hunds_case_a_state(label=self.label, Lambda=self.Lambda, 
                                                    S=self.S, Sigma=Sigma, J=self.J, 
                                                    Omega=Sigma+self.Lambda, I=self.I, F=self.F, mF=self.mF)
                    self.Hunds_case_a_list.append(d)


class SrFacStarkShift:
    def __init__(self, state, LaserWavelength_nm=1064, LaserIntensity_kW_invcm2=0, LaserPolarization=0, print_dipole_moment=False, print_polarizability=True, print_stark_shift=True, print_scattering_rate=True):
        self.state = state
        self.LaserWavelength_nm = LaserWavelength_nm
        self.LaserIntensity_kW_invcm2 = LaserIntensity_kW_invcm2

        self.define_constants()
        self.calculate_dipole_moment(print_result=print_dipole_moment)
        self.calculate_polarizability(print_result=print_polarizability)
        self.calculate_stark_shift(print_result=print_stark_shift)
        self.calculate_scattering_rate(LaserPolarization=LaserPolarization, print_result=print_scattering_rate)

    def define_constants(self):
        # spectroscopic data are from Thomas' Matlab code
        # define electronic state energies in wavenumber (1/wavelength, cm^{-1})
        XSigmaEnergy_invcm = 0
        APiOneHalfEnergy_invcm = 15075.6122
        APiThreeHalvesEnergy_invcm = 15357.0736
        BSigmaEnergy_invcm = 17267.4465
        CPiOneHalfEnergy_invcm = 27384.67 - 57.9048/2
        CPiThreeHalvesEnergy_invcm = 27384.67 + 57.9048/2
        DSigmaEnergy_invcm = 27773.8
        FSigmaEnergy_invcm = 32823.5
        GPiEnergy_invcm = 34808.9275
        HEnergy_invcm = 42345

        # define electronic state lifetimes in ns
        self.APiOneHalfLifetime_ns = 24.1
        APiThreeHalvesLifetime_ns = 22.6
        BSigmaLifetime_ns = 25.5
        CPiOneHalfLifetime_ns = 66.3
        DSigmaLifetime_ns = 219

        # define fundamental constants in atomic units
        self.FineStructureConstant = 1/137.035999084
        self.SpeedOfLight = 1/self.FineStructureConstant
        self.e = 1 # electron charge
        self.hbar = 1
        self.DebyeToAtomicUnit = 0.3934303 # Debye in atomic units
        nmToAtomicUnit = 1/0.0529177210903
        invcmToAtomicUnit = 1e-7/nmToAtomicUnit*self.SpeedOfLight*self.hbar*2*np.pi
        nsToAtomicUnit = 1/(2.4188843265857e-8)

        LaserWavelength = self.LaserWavelength_nm*nmToAtomicUnit
        LaserEnergy = self.SpeedOfLight/LaserWavelength*self.hbar*2*np.pi

        # define electronic state energies in atomic units (Hartree energy)
        self.XSigmaEnergy = XSigmaEnergy_invcm*invcmToAtomicUnit
        self.APiOneHalfEnergy = APiOneHalfEnergy_invcm*invcmToAtomicUnit
        self.APiThreeHalvesEnergy = APiThreeHalvesEnergy_invcm*invcmToAtomicUnit
        self.BSigmaEnergy = BSigmaEnergy_invcm*invcmToAtomicUnit
        self.CPiOneHalfEnergy = CPiOneHalfEnergy_invcm*invcmToAtomicUnit
        self.CPiThreeHalvesEnergy = CPiThreeHalvesEnergy_invcm*invcmToAtomicUnit
        self.DSigmaEnergy = DSigmaEnergy_invcm*invcmToAtomicUnit
        self.FSigmaEnergy = FSigmaEnergy_invcm*invcmToAtomicUnit
        self.GPiEnergy = GPiEnergy_invcm*invcmToAtomicUnit
        self.HEnergy = HEnergy_invcm*invcmToAtomicUnit

        # define detunings relative to each electronic state, co-rotating term in atomic units
        self.DetuningAPiOneHalf = self.APiOneHalfEnergy - LaserEnergy
        self.DetuningAPiThreeHalves = self.APiThreeHalvesEnergy - LaserEnergy
        self.DetuningBSigma = self.BSigmaEnergy - LaserEnergy
        self.DetuningCPiOneHalf = self.CPiOneHalfEnergy - LaserEnergy
        self.DetuningCPiThreeHalves = self.CPiThreeHalvesEnergy - LaserEnergy
        self.DetuningDSigma = self.DSigmaEnergy - LaserEnergy
        self.DetuningFSigma = self.FSigmaEnergy - LaserEnergy
        self.DetuningGPi = self.GPiEnergy - LaserEnergy
        self.DetuningH = self.HEnergy - LaserEnergy

        # define detunings relative to each electronic state, counter-rotating term in atomic units
        self.DetuningAPiOneHalfCounterRotate = self.APiOneHalfEnergy + LaserEnergy
        self.DetuningAPiThreeHalvesCounterRotate = self.APiThreeHalvesEnergy + LaserEnergy
        self.DetuningBSigmaCounterRotate = self.BSigmaEnergy + LaserEnergy
        self.DetuningCPiOneHalfCounterRotate = self.CPiOneHalfEnergy + LaserEnergy
        self.DetuningCPiThreeHalvesCounterRotate = self.CPiThreeHalvesEnergy + LaserEnergy
        self.DetuningDSigmaCounterRotate = self.DSigmaEnergy + LaserEnergy
        self.DetuningFSigmaCounterRotate = self.FSigmaEnergy + LaserEnergy
        self.DetuningGPiCounterRotate = self.GPiEnergy + LaserEnergy
        self.DetuningHCounterRotate = self.HEnergy + LaserEnergy

        # define electronic state lifetimes in atomic units
        self.APiOneHalfLifetime = self.APiOneHalfLifetime_ns*nsToAtomicUnit
        self.APiThreeHalvesLifetime = APiThreeHalvesLifetime_ns*nsToAtomicUnit
        self.BSigmaLifetime = BSigmaLifetime_ns*nsToAtomicUnit
        self.CPiOneHalfLifetime = CPiOneHalfLifetime_ns*nsToAtomicUnit
        self.DSigmaLifetime = DSigmaLifetime_ns*nsToAtomicUnit

    def calculate_dipole_moment(self, print_result=True):
        self.excited_state_list = []

        # A Pi one half state
        TransitionFreq = (self.APiOneHalfEnergy - self.XSigmaEnergy)/self.hbar
        APiOneHalfDipole = np.sqrt(3*(self.SpeedOfLight**2)*(self.e**2)/4/self.FineStructureConstant/(TransitionFreq**3)/self.APiOneHalfLifetime)
        self.excited_state_list.append({"label":"APiOneHalf", "Lambda":1, "Sigma":-1/2, "transition dipole":APiOneHalfDipole, 
                                        "detuning co-rotate":self.DetuningAPiOneHalf, "detuning counter-rotate":self.DetuningAPiOneHalfCounterRotate})
        self.excited_state_list.append({"label":"APiOneHalf", "Lambda":-1, "Sigma":1/2, "transition dipole":APiOneHalfDipole,
                                        "detuning co-rotate":self.DetuningAPiOneHalf, "detuning counter-rotate":self.DetuningAPiOneHalfCounterRotate})
        
        # A Pi three halves state
        TransitionFreq = (self.APiThreeHalvesEnergy - self.XSigmaEnergy)/self.hbar
        APiThreeHalvesDipole = np.sqrt(3*(self.SpeedOfLight**2)*(self.e**2)/4/self.FineStructureConstant/(TransitionFreq**3)/self.APiThreeHalvesLifetime)
        self.excited_state_list.append({"label":"APiThreeHalves", "Lambda":1, "Sigma":1/2, "transition dipole":APiThreeHalvesDipole, 
                                        "detuning co-rotate":self.DetuningAPiThreeHalves, "detuning counter-rotate":self.DetuningAPiThreeHalvesCounterRotate})
        self.excited_state_list.append({"label":"APiThreeHalves", "Lambda":-1, "Sigma":-1/2, "transition dipole":APiThreeHalvesDipole, 
                                        "detuning co-rotate":self.DetuningAPiThreeHalves, "detuning counter-rotate":self.DetuningAPiThreeHalvesCounterRotate})

        # B Sigma state
        TransitionFreq = (self.BSigmaEnergy - self.XSigmaEnergy)/self.hbar
        BSigmaDipole = np.sqrt(3*(self.SpeedOfLight**2)*(self.e**2)/4/self.FineStructureConstant/(TransitionFreq**3)/self.BSigmaLifetime)
        self.excited_state_list.append({"label":"BSigma", "Lambda":0, "Sigma":1/2, "transition dipole":BSigmaDipole, 
                                        "detuning co-rotate":self.DetuningBSigma, "detuning counter-rotate":self.DetuningBSigmaCounterRotate})
        self.excited_state_list.append({"label":"BSigma", "Lambda":0, "Sigma":-1/2, "transition dipole":BSigmaDipole, 
                                        "detuning co-rotate":self.DetuningBSigma, "detuning counter-rotate":self.DetuningBSigmaCounterRotate})
        
        # C Pi state
        TransitionFreq = (self.CPiOneHalfEnergy - self.XSigmaEnergy)/self.hbar
        CPiOneHalfDipole = np.sqrt(3*(self.SpeedOfLight**2)*(self.e**2)/4/self.FineStructureConstant/(TransitionFreq**3)/self.CPiOneHalfLifetime)
        self.excited_state_list.append({"label":"CPi", "Lambda":1, "Sigma":-1/2, "transition dipole":CPiOneHalfDipole, 
                                        "detuning co-rotate":self.DetuningCPiOneHalf, "detuning counter-rotate":self.DetuningCPiOneHalfCounterRotate})
        self.excited_state_list.append({"label":"CPi", "Lambda":-1, "Sigma":1/2, "transition dipole":CPiOneHalfDipole, 
                                        "detuning co-rotate":self.DetuningCPiOneHalf, "detuning counter-rotate":self.DetuningCPiOneHalfCounterRotate})

        # assume CPi three halves states have the same dipole moment as CPi one half states
        self.excited_state_list.append({"label":"CPi", "Lambda":1, "Sigma":1/2, "transition dipole":CPiOneHalfDipole, 
                                        "detuning co-rotate":self.DetuningCPiThreeHalves, "detuning counter-rotate":self.DetuningCPiThreeHalvesCounterRotate}) 
        self.excited_state_list.append({"label":"CPi", "Lambda":-1, "Sigma":-1/2, "transition dipole":CPiOneHalfDipole,
                                        "detuning co-rotate":self.DetuningCPiThreeHalves, "detuning counter-rotate":self.DetuningCPiThreeHalvesCounterRotate})

        # D Sigma state
        TransitionFreq = (self.DSigmaEnergy - self.XSigmaEnergy)/self.hbar
        DSigmaDipole = np.sqrt(3*(self.SpeedOfLight**2)*(self.e**2)/4/self.FineStructureConstant/(TransitionFreq**3)/self.DSigmaLifetime)
        self.excited_state_list.append({"label":"DSigma", "Lambda":0, "Sigma":1/2, "transition dipole":DSigmaDipole, 
                                        "detuning co-rotate":self.DetuningDSigma, "detuning counter-rotate":self.DetuningDSigmaCounterRotate})
        self.excited_state_list.append({"label":"DSigma", "Lambda":0, "Sigma":-1/2, "transition dipole":DSigmaDipole, 
                                        "detuning co-rotate":self.DetuningDSigma, "detuning counter-rotate":self.DetuningDSigmaCounterRotate})

        # F Sigma state
        FSigmaDipole = 0.4
        self.excited_state_list.append({"label":"FSigma", "Lambda":0, "Sigma":1/2, "transition dipole":FSigmaDipole, 
                                        "detuning co-rotate":self.DetuningFSigma, "detuning counter-rotate":self.DetuningFSigmaCounterRotate})
        self.excited_state_list.append({"label":"FSigma", "Lambda":0, "Sigma":-1/2, "transition dipole":FSigmaDipole, 
                                        "detuning co-rotate":self.DetuningFSigma, "detuning counter-rotate":self.DetuningFSigmaCounterRotate})

        # G Pi state
        GPiDipole = 0.5
        self.excited_state_list.append({"label":"GPi", "Lambda":1, "Sigma":-1/2, "transition dipole":GPiDipole, 
                                        "detuning co-rotate":self.DetuningGPi, "detuning counter-rotate":self.DetuningGPiCounterRotate})
        self.excited_state_list.append({"label":"GPi", "Lambda":-1, "Sigma":1/2, "transition dipole":GPiDipole, 
                                        "detuning co-rotate":self.DetuningGPi, "detuning counter-rotate":self.DetuningGPiCounterRotate})
        self.excited_state_list.append({"label":"GPi", "Lambda":1, "Sigma":1/2, "transition dipole":GPiDipole, 
                                        "detuning co-rotate":self.DetuningGPi, "detuning counter-rotate":self.DetuningGPiCounterRotate})
        self.excited_state_list.append({"label":"GPi", "Lambda":-1, "Sigma":-1/2, "transition dipole":GPiDipole, 
                                        "detuning co-rotate":self.DetuningGPi, "detuning counter-rotate":self.DetuningGPiCounterRotate})

        # H state
        HDipole = 2.2934*self.DebyeToAtomicUnit

        if print_result:
            print("Transition dipole moments from X state.")
            print("---------------------------------------")
            print("A Pi one half state: {:.3f} a.u.".format(APiOneHalfDipole))
            print("A Pi three halves state: {:.3f} a.u.".format(APiThreeHalvesDipole))
            print("B Sigma state: {:.3f} a.u.".format(BSigmaDipole))
            print("C Pi one half state: {:.3f} a.u.".format(CPiOneHalfDipole))
            print("D Sigma state: {:.3f} a.u.".format(DSigmaDipole))
            print("F Sigma state: {:.3f} a.u.".format(FSigmaDipole))
            print("G Pi state: {:.3f} a.u.".format(GPiDipole))
            print("H state: {:.3f} a.u.".format(HDipole))
            print("")

    def calculate_polarizability(self, print_result=True):
        self.polarizability_dict = {}
        for K in [0, 1, 2]:
            partial_polarizability = 0
            for excited_state in self.excited_state_list:
                # if excited_state["label"] not in ["APiOneHalf", "APiThreeHalves"]:
                #     continue

                Lambda_prime = excited_state["Lambda"]
                Sigma_prime = excited_state["Sigma"]
                Omega_prime = Lambda_prime + Sigma_prime

                for J_prime in np.arange(np.abs(Omega_prime), self.state.J+10):
                    # just pick a big number for upper limit of J_prime

                    for F_prime in np.arange(np.abs(J_prime-self.state.I), J_prime+self.state.I+1):
                        partial_dipole = 0
                        for d in self.state.Hunds_case_a_list:
                            state_prime = Hunds_case_a_state(label=excited_state["label"], Lambda=Lambda_prime, S=self.state.S, 
                                                            Sigma=Sigma_prime, J=J_prime, Omega=Omega_prime, I=self.state.I, F=F_prime, mF=0) # mF doesn't matter her, just put a random number
                            partial_dipole += d["coefficient"]*(d["state"].calculate_dipole_moment(dipole_moment_Lambda=excited_state["transition dipole"], state_prime=state_prime))
                        partial_polarizability += (((-1)**(self.state.F+F_prime))*float(wigner_6j(1, 1, K, self.state.F, self.state.F, F_prime))*(partial_dipole**2)
                                                    *(1/excited_state["detuning co-rotate"]+((-1)**K)/excited_state["detuning counter-rotate"]))
            
            partial_polarizability *= np.sqrt(2*K+1)*(-1)**K
            if K == 0:
                polarizability = (-1)*partial_polarizability/np.sqrt(3*(2*self.state.F+1))
                self.polarizability_dict["scalar polarizability"] = polarizability
            elif K == 1:
                polarizability = partial_polarizability*np.sqrt(2*self.state.F/(self.state.F+1)/(2*self.state.F+1))
                self.polarizability_dict["vector polarizability"] = polarizability
            elif K == 2:
                polarizability = partial_polarizability*np.sqrt(2*self.state.F*(2*self.state.F-1)/3/(self.state.F+1)/(2*self.state.F+1)/(2*self.state.F+3))
                self.polarizability_dict["tensor polarizability"] = polarizability
            else:
                print(f"Incorrect rank K: {K}")
                return

        if print_result:
            print(f"Polarizability of state {self.state.label} (N={self.state.N}, J={self.state.J}, F={self.state.F}) at {self.LaserWavelength_nm} nm.")
            print("----------------------------------------------------------")
            for pol_type, pol in self.polarizability_dict.items():
                print(pol_type+": {:.2f} a.u.".format(pol))
            print("")

    def calculate_stark_shift(self, print_result=True):
        permittivity = 1/(4*np.pi) # permittivity in atomic units
        kw_invcm2_ToAtomicUnit = 1/((4.3597447222071e-18)**2/(1.054571817e-34)/1e3/(5.29177210903e-9)**2)
        HartreeTouK = 315775.02480407e6 # convert Hartree to uK
        HartreeToMHz =  6579.683920502e6 # convert Hartree to MHz
        LaserIntensity = self.LaserIntensity_kW_invcm2*kw_invcm2_ToAtomicUnit # Laser intensity in atomic units
        prefactor = LaserIntensity*2/permittivity/self.SpeedOfLight/4 # prefactor for converting polarizability to ac Stark shift

        ScalarShift = prefactor*self.polarizability_dict["scalar polarizability"]
        ScalarShift_uK = ScalarShift*HartreeTouK
        ScalarShift_MHz = ScalarShift*HartreeToMHz

        self.StarkShift_dict = {}
        self.StarkShift_dict["scalar shift uK"] = ScalarShift_uK
        self.StarkShift_dict["scalar shift MHz"] = ScalarShift_MHz

        if print_result:
            print(f"AC Stark Shift of state {self.state.label} (N={self.state.N}, J={self.state.J}, F={self.state.F}) at {self.LaserWavelength_nm} nm.")
            print("----------------------------------------------------------")
            print("scalar shift: {:.2f} uK ({:.2f} MHz)".format(ScalarShift_uK, ScalarShift_MHz))
            print("")

    def calculate_scattering_rate(self, LaserPolarization=0, print_result=True):
        # HartreeTouK = 315775.02480407e6 # convert Hartree to uK
        # u = self.StarkShift_dict["scalar shift uK"]/HartreeTouK
        # Delta = self.DetuningAPiOneHalf
        # Gamma_sc = u/Delta*(1/(self.APiOneHalfLifetime_ns/1e9)) # scattering rate in 1/s
        # self.scattering_rate_invs = np.abs(Gamma_sc)

        # hbar = 1.05457182e-34 # SI units
        # k = 2*np.pi/(self.LaserWavelength_nm/1e9) # 1/m
        # m = 106.62*1.66053906660e-27 # kg
        # kB = 1.380649e-23 # J/K
        # T_rec = ((hbar*k)**2)/m/kB # recoil temperature, K
        # self.heating_rate_nK_s = T_rec*self.scattering_rate_invs/3*1e9 # nK/s

        LaserEnergy = (1/(self.LaserWavelength_nm/1e9)*299792458*2*np.pi)/(4.134137336e16) # laser energy in atomic units
        permittivity = 1/(4*np.pi) # permittivity in atomic units
        kw_invcm2_ToAtomicUnit = 1/((4.3597447222071e-18)**2/(1.054571817e-34)/1e3/(5.29177210903e-9)**2)
        AtomicUnit_To_seconds = 2.4188843265857e-17 # seconds
        LaserIntensity = self.LaserIntensity_kW_invcm2*kw_invcm2_ToAtomicUnit # Laser intensity in atomic units
        prefactor = LaserIntensity*LaserEnergy**3/(6*np.pi)/permittivity**2/self.SpeedOfLight**4 # prefactor for converting polarizability to ac Stark shift

        initial_state = self.state # assume the initial state to be self.state
        final_state = self.state # assume the final state is in ground electronic state, we only use J from self.state, not F
        Jprime = final_state.J
        I = initial_state.I
        J = initial_state.J
        F = initial_state.F
        mF = initial_state.mF

        s = 0



class Rb_coupled_state:
    def __init__(self, label='', energy=0, J=0, reduced_matrix_element=0):
        self.label = label
        self.energy = energy # in 1/cm
        self.energy_au = (energy*1e2*299792458*2*np.pi)/(4.134137336e16) # energy in atomic units, 1 a.u = 4.134137336*10^16 Hz (angular frequency)
        self.J = J
        self.reduced_matrix_element = reduced_matrix_element # in atomic unit, namely ea_0, where a_0 is Bohr radius  


# adapted from my august 2021 code, in slighly different format from SrF calculation
class Rb_ground_state:
    def __init__(self, label='5S1_2 state', energy=0, J=1/2, I=3/2, F=2, mF=0):
        self.label = label
        self.energy = energy # in 1/cm
        self.energy_au = (energy*1e2*299792458*2*np.pi)/(4.134137336e16) # energy in atomic units, 1 a.u = 4.134137336*10^16 Hz (angular frequency)
        self.J = J
        self.I = I
        self.F = F
        self.mF = mF

        self.coupled_state_list = []
        self.coupled_state_list.append(Rb_coupled_state(label='5p1_2', energy=12578.950, J=1/2, reduced_matrix_element=4.227))
        self.coupled_state_list.append(Rb_coupled_state(label='5p3_2', energy=12816.545, J=3/2, reduced_matrix_element=5.977))
        self.coupled_state_list.append(Rb_coupled_state(label='6p1_2', energy=23715.081, J=1/2, reduced_matrix_element=0.342))
        self.coupled_state_list.append(Rb_coupled_state(label='6p3_2', energy=23792.591, J=3/2, reduced_matrix_element=0.553))
        self.coupled_state_list.append(Rb_coupled_state(label='7p1_2', energy=27835.02, J=1/2, reduced_matrix_element=0.118))
        self.coupled_state_list.append(Rb_coupled_state(label='7p3_2', energy=27870.11, J=3/2, reduced_matrix_element=0.207))
        self.coupled_state_list.append(Rb_coupled_state(label='8p1_2', energy=29834.94, J=1/2, reduced_matrix_element=0.061))
        self.coupled_state_list.append(Rb_coupled_state(label='8p3_2', energy=29853.79, J=3/2, reduced_matrix_element=0.114))
        self.coupled_state_list.append(Rb_coupled_state(label='9p1_2', energy=30958.91, J=1/2, reduced_matrix_element=0.046))
        self.coupled_state_list.append(Rb_coupled_state(label='9p3_2', energy=30970.19, J=3/2, reduced_matrix_element=0.074))

        # core polarizability at 1064nm, in atomic units. Wavelength dependence is weak. Only contributes to scalar polarizability. The same for all Rb states 
        self.polarizability_core = 9.3

class RbacStarkShift:
    def __init__(self, state, LaserWavelength_nm=1064, LaserIntensity_kW_invcm2=0, LaserPolarization=0,
                 print_polarizability=True, print_stark_shift=True, print_scattering_rate=True):
        self.state = state
        self.LaserWavelength_nm = LaserWavelength_nm
        self.LaserIntensity_kW_invcm2 = LaserIntensity_kW_invcm2

        self.calculate_polarizability(print_result=print_polarizability)
        self.calculate_stark_shift(print_result=print_stark_shift)
        self.calculate_scattering_rate(LaserPolarization=LaserPolarization, print_result=print_scattering_rate)

    def calculate_polarizability(self, print_result=True):
        LaserEnergy = (1/(self.LaserWavelength_nm/1e9)*299792458*2*np.pi)/(4.134137336e16) # laser energy in atomic units
        J = self.state.J
        I = self.state.I
        F = self.state.F

        self.polarizability_dict = {}

        for K in [0, 1, 2]:
            pol_K_J = 0
            for state in self.state.coupled_state_list:
                pol_K_J += (-1)**(state.J+J)*wigner_6j(1, 1, K, J, J, state.J)*(state.reduced_matrix_element**2)*(1/(state.energy_au-self.state.energy_au-LaserEnergy)+
                                                                                                            (-1)**K/(state.energy_au-self.state.energy_au+LaserEnergy))
            pol_K_J *= np.sqrt(2*K+1)*(-1)**K
            pol_K_F = (-1)**(F+J+K+I)*(2*F+1)*wigner_6j(J, F, I, F, J, K)*pol_K_J

            if K == 0:
                pol_F = (-1/np.sqrt(3*(2*F+1)))*pol_K_F
                pol_F += self.state.polarizability_core
                self.polarizability_dict["scalar polarizability"] = float(pol_F)
            elif K == 1:
                pol_F = np.sqrt(2*F/(F+1)/(2*F+1))*pol_K_F
                self.polarizability_dict["vector polarizability"] = float(pol_F)
            elif K == 2:
                pol_F = np.sqrt(2*F*(2*F-1)/3/(F+1)/(2*F+1)/(2*F+3))*pol_K_F
                self.polarizability_dict["tensor polarizability"] = float(pol_F)
            else:
                print("Invalid rank.")
                return

        if print_result:
            print(f"Polarizability of state {self.state.label} (J={self.state.J}, F={self.state.F}) at {self.LaserWavelength_nm} nm.")
            print("-----------------------------------------------------")
            for pol_type, pol in self.polarizability_dict.items():
                print(pol_type+": {:.2f} a.u.".format(pol))
            print("")

    def calculate_stark_shift(self, print_result=True):
        permittivity = 1/(4*np.pi) # permittivity in atomic units
        kw_invcm2_ToAtomicUnit = 1/((4.3597447222071e-18)**2/(1.054571817e-34)/1e3/(5.29177210903e-9)**2)
        HartreeTouK = 315775.02480407e6 # convert Hartree to uK
        HartreeToMHz =  6579.683920502e6 # convert Hartree to MHz
        SpeedOfLight = 137.035999084 # in atomic units
        LaserIntensity = self.LaserIntensity_kW_invcm2*kw_invcm2_ToAtomicUnit # Laser intensity in atomic units
        prefactor = LaserIntensity*2/permittivity/SpeedOfLight/4 # prefactor for converting polarizability to ac Stark shift

        ScalarShift = prefactor*self.polarizability_dict["scalar polarizability"]
        ScalarShift_uK = ScalarShift*HartreeTouK
        ScalarShift_MHz = ScalarShift*HartreeToMHz

        self.StarkShift_dict = {}
        self.StarkShift_dict["scalar shift uK"] = ScalarShift_uK
        self.StarkShift_dict["scalar shift MHz"] = ScalarShift_MHz

        if print_result:
            print(f"AC Stark Shift of state {self.state.label} (J={self.state.J}, F={self.state.F}) at {self.LaserWavelength_nm} nm.")
            print("----------------------------------------------------------")
            print("scalar shift: {:.2f} uK ({:.2f} MHz)".format(ScalarShift_uK, ScalarShift_MHz))
            print("")

    def calculate_scattering_rate(self, LaserPolarization=0, print_result=True):     
        LaserEnergy = (1/(self.LaserWavelength_nm/1e9)*299792458*2*np.pi)/(4.134137336e16) # laser energy in atomic units
        permittivity = 1/(4*np.pi) # permittivity in atomic units
        kw_invcm2_ToAtomicUnit = 1/((4.3597447222071e-18)**2/(1.054571817e-34)/1e3/(5.29177210903e-9)**2)
        AtomicUnit_To_seconds = 2.4188843265857e-17 # seconds
        SpeedOfLight = 137.035999084 # in atomic units
        LaserIntensity = self.LaserIntensity_kW_invcm2*kw_invcm2_ToAtomicUnit # Laser intensity in atomic units
        prefactor = LaserIntensity*LaserEnergy**3/(6*np.pi)/permittivity**2/SpeedOfLight**4 # prefactor for converting polarizability to ac Stark shift

        initial_state = self.state # assume the initial state to be self.state
        final_state = self.state # assume the final state is in ground electronic state, we only use J from self.state, not F
        Jprime = final_state.J
        I = initial_state.I
        J = initial_state.J
        F = initial_state.F
        mF = initial_state.mF

        s = 0
        for Fprime in np.arange(np.abs(Jprime-I), Jprime+I+1):
            for mFprime in np.arange(-Fprime, Fprime+1):
                for p in [-1, 0, 1]:
                    a = 0

                    for intermediate_state in self.state.coupled_state_list:
                        Jdoubleprime = intermediate_state.J

                        for Fdoubleprime in np.arange(np.abs(Jdoubleprime-I), Jdoubleprime+I+1):
                            for mFdoubleprime in np.arange(-Fdoubleprime, Fdoubleprime+1):
                                b = 1

                                b *= wigner_eckart_coefficient(Fprime, mFprime, 1, p, Fdoubleprime, mFdoubleprime)
                                b *= spectator_theorem_coefcient(Jprime, I, Fprime, 1, Jdoubleprime, I, Fdoubleprime)

                                b *= wigner_eckart_coefficient(Fdoubleprime, mFdoubleprime, 1, LaserPolarization, F, mF)
                                b *= spectator_theorem_coefcient(Jdoubleprime, I, Fdoubleprime, 1, J, I, F)

                                b *= (-1)**(Jdoubleprime-J)*intermediate_state.reduced_matrix_element**2

                                b /= (intermediate_state.energy_au-LaserEnergy)

                                a += b

                                b = 1

                                b *= wigner_eckart_coefficient(Fprime, mFprime, 1, LaserPolarization, Fdoubleprime, mFdoubleprime)
                                b *= spectator_theorem_coefcient(Jprime, I, Fprime, 1, Jdoubleprime, I, Fdoubleprime)

                                b *= wigner_eckart_coefficient(Fdoubleprime, mFdoubleprime, 1, p, F, mF)
                                b *= spectator_theorem_coefcient(Jdoubleprime, I, Fdoubleprime, 1, J, I, F)

                                b *= (-1)**(Jdoubleprime-J)*intermediate_state.reduced_matrix_element**2

                                b /= (intermediate_state.energy_au+LaserEnergy)

                                a += b

                    s += np.abs(a)**2

        scattering_rate_au = prefactor*s
        scattering_rate = scattering_rate_au/AtomicUnit_To_seconds
        self.scattering_rate_invs = scattering_rate

        hbar = 1.05457182e-34 # SI units
        k = 2*np.pi/(self.LaserWavelength_nm/1e9) # 1/m
        m = 86.9*1.66053906660e-27 # kg
        kB = 1.380649e-23 # J/K
        T_rec = ((hbar*k)**2)/m/kB # recoil temperature, K
        self.heating_rate_nK_s = T_rec*self.scattering_rate_invs/3*1e9 # nK/s

        if print_result:
            print(f"Photon scattering rate of state {self.state.label} (J={self.state.J}, F={self.state.F}, mF={mF}) in {self.LaserWavelength_nm} nm (polarization={LaserPolarization}) ODT.")
            print("----------------------------------------------------------")
            print("Scattering rate: {:.2f} 1/s ".format(self.scattering_rate_invs))
            print("Heatinging rate: {:.0f} nK/s ".format(self.heating_rate_nK_s))
            print("")
