import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import os, time
from classes import *

def trap_gravity(x, optical_trap_depth, sigma, m):
    """
    calculate ODT trap potential, as a function of position x, under the influence of gravity
    """

    # x is position, in um
    # optical_trap_depth in uK
    # sigma is sigma of ODT gaussian distribution, in um
    # m is mass in atomic mass unit

    kB = 1.380649e-23 # SI units
    x *= 1.e-6 # convert to m
    sigma *= 1.e-6 # convert to m
    m *= 1.66053906660e-27 # kg
    g = 9.8 # kg*m/s^2

    r = -(optical_trap_depth/1.e6)*kB*np.exp(-x**2/2/sigma**2) + m*g*x

    return r/kB*1e6 # in uK

def trap_gravity_derivative(x, optical_trap_depth, sigma, m):
    """
    calculate the position derivative of ODT trap potential, under the influence of gravity
    """

    # x is position, in um
    # optical_trap_depth in uK
    # sigma is sigma of ODT gaussian distribution, in um
    # m is mass in atomic mass unit

    kB = 1.380649e-23 # SI units
    x = x*1.e-6 # convert to m, x is an array
    sigma *= 1.e-6 # convert to m
    m *= 1.66053906660e-27 # kg
    g = 9.8 # kg*m/s^2

    r = (optical_trap_depth/1.e6)*kB*np.exp(-x[0]**2/2/sigma**2)*(x[0]/sigma**2) + m*g

    return r/kB # in K/m


def BichromaticTrap(target_shift_Rb, target_shift_SrF, wavelength_laser2, radius, print_result=True):
    SrF_mass = 106.62
    Rb_mass = 86.9

    # make initial guess and calculate trap depth
    # ----------------------------------------------------------------------------------------------------------------
    power_1064_initial = 1 # laser power in W
    intensity_1064 = (power_1064_initial/2/np.pi/(radius/2)**2)/1e3*1e8 # peak Laser intensity in kW/cm^2

    power_laser2_initial = 1 # laser power in W
    intensity_laser2 = (power_laser2_initial/2/np.pi/(radius/2)**2)/1e3*1e8 # peak Laser intensity in kW/cm^2

    SrF_state = Hunds_case_b_state(label="SrF XSigma", Lambda=0, N=0, S=1/2, J=1/2, I=1/2, F=0, mF=0)
    Rb_state = Rb_ground_state(label="Rb 5S", energy=0, J=1/2, I=3/2, F=1, mF=0)

    SrF_shift_1064_initial = SrFacStarkShift(SrF_state, LaserWavelength_nm=1064, LaserIntensity_kW_invcm2=intensity_1064, print_dipole_moment=False, print_polarizability=False, print_stark_shift=False, print_scattering_rate=False)
    Rb_shift_1064_initial = RbacStarkShift(Rb_state, LaserWavelength_nm=1064, LaserIntensity_kW_invcm2=intensity_1064, print_polarizability=False, print_stark_shift=False, print_scattering_rate=False)

    SrF_shift_laser2_initial = SrFacStarkShift(SrF_state, LaserWavelength_nm=wavelength_laser2, LaserIntensity_kW_invcm2=intensity_laser2, print_dipole_moment=False, print_polarizability=False, print_stark_shift=False, print_scattering_rate=False) 
    Rb_shift_laser2_initial = RbacStarkShift(Rb_state, LaserWavelength_nm=wavelength_laser2, LaserIntensity_kW_invcm2=intensity_laser2, print_polarizability=False, print_stark_shift=False, print_scattering_rate=False)
    # ------------------------------------------------------------------------------------------------------------------------


    # calculate required power to get target depth
    #---------------------------------------------------------------------------------------------------------------------------
    a = np.array([[Rb_shift_1064_initial.StarkShift_dict["scalar shift uK"]/power_1064_initial, Rb_shift_laser2_initial.StarkShift_dict["scalar shift uK"]/power_laser2_initial], 
                    [SrF_shift_1064_initial.StarkShift_dict["scalar shift uK"]/power_1064_initial, SrF_shift_laser2_initial.StarkShift_dict["scalar shift uK"]/power_laser2_initial]])
    b = np.array([target_shift_Rb, target_shift_SrF])
    power_1064_final, power_laser2_final = np.linalg.solve(a, b)
    #---------------------------------------------------------------------------------------------------------------------------


    # re-calculate stark shift for two species to confirm they have the same trap depth, and reach the target depth
    #--------------------------------------------------------------------------------------------------------------------------
    SrF_shift_1064_final = SrFacStarkShift(SrF_state, LaserWavelength_nm=1064, LaserIntensity_kW_invcm2=intensity_1064*(power_1064_final/power_1064_initial), print_dipole_moment=False, print_polarizability=False, print_stark_shift=False, print_scattering_rate=False)
    Rb_shift_1064_final = RbacStarkShift(Rb_state, LaserWavelength_nm=1064, LaserIntensity_kW_invcm2=intensity_1064*(power_1064_final/power_1064_initial), print_polarizability=False, print_stark_shift=False, print_scattering_rate=False)

    SrF_shift_laser2_final = SrFacStarkShift(SrF_state, LaserWavelength_nm=wavelength_laser2, LaserIntensity_kW_invcm2=intensity_laser2*(power_laser2_final/power_laser2_initial), print_dipole_moment=False, print_polarizability=False, print_stark_shift=False, print_scattering_rate=False)
    Rb_shift_laser2_final = RbacStarkShift(Rb_state, LaserWavelength_nm=wavelength_laser2, LaserIntensity_kW_invcm2=intensity_laser2*(power_laser2_final/power_laser2_initial), print_polarizability=False, print_stark_shift=False, print_scattering_rate=False)

    SrF_shift_final = SrF_shift_1064_final.StarkShift_dict["scalar shift uK"] + SrF_shift_laser2_final.StarkShift_dict["scalar shift uK"]
    SrF_scattering_rate = SrF_shift_1064_final.scattering_rate_invs + SrF_shift_laser2_final.scattering_rate_invs
    SrF_heating_rate = SrF_shift_1064_final.heating_rate_nK_s + SrF_shift_laser2_final.heating_rate_nK_s

    Rb_shift_final = Rb_shift_1064_final.StarkShift_dict["scalar shift uK"] + Rb_shift_laser2_final.StarkShift_dict["scalar shift uK"]
    Rb_scattering_rate = Rb_shift_1064_final.scattering_rate_invs + Rb_shift_laser2_final.scattering_rate_invs
    Rb_heating_rate = Rb_shift_1064_final.heating_rate_nK_s + Rb_shift_laser2_final.heating_rate_nK_s
    #--------------------------------------------------------------------------------------------------------------------------


    # calculte trap depth reduce caused by gravity
    # ------------------------------------------------------------------------------------------------------------------------
    f = lambda x: trap_gravity_derivative(x, optical_trap_depth=SrF_shift_final, sigma=radius/2, m=SrF_mass)
    root1 = fsolve(f, [0])[0] # um
    root2 = fsolve(f, [-radius])[0] # um
    SrF_trap_depth = trap_gravity(root2, optical_trap_depth=SrF_shift_final, sigma=radius/2, m=SrF_mass) - trap_gravity(root1, optical_trap_depth=SrF_shift_final, sigma=radius/2, m=SrF_mass)

    f = lambda x: trap_gravity_derivative(x, optical_trap_depth=Rb_shift_final, sigma=radius/2, m=Rb_mass)
    root1 = fsolve(f, [0])[0] # um
    root2 = fsolve(f, [-radius])[0] # um
    Rb_trap_depth = trap_gravity(root2, optical_trap_depth=Rb_shift_final, sigma=radius/2, m=Rb_mass) - trap_gravity(root1, optical_trap_depth=Rb_shift_final, sigma=radius/2, m=Rb_mass)
    # ------------------------------------------------------------------------------------------------------------------------


    if np.abs(SrF_shift_final - target_shift_SrF) > 1e-6:
        print('')
        print("SrF Stark shift doesn't reach the target value.")
        print("SrF target shift: {:.2f} uK".format(target_shift_SrF))
        print("SrF actual shift: {:.2f} uK".format(SrF_shift_final))
        return

    if np.abs(Rb_shift_final - target_shift_Rb) > 1e-6:
        print('')
        print("Rb Stark shift doesn't reach the target value.")
        print("Rb target shift: {:.2f} uK".format(target_shift_Rb))
        print("Rb actual shift: {:.2f} uK".format(Rb_shift_final))
        return

    if print_result:
        print('')
        print("SrF Stark shift: {:.2f} uk".format(SrF_shift_final))
        print("SrF Stark shift from 1064 nm laser: {:.2f} uk".format(SrF_shift_1064_final.StarkShift_dict["scalar shift uK"]))
        print("SrF Stark shift from laser2: {:.2f} uk".format(SrF_shift_laser2_final.StarkShift_dict["scalar shift uK"]))
        print("SrF trap depth w/ gravity: {:.2f} uk".format(SrF_trap_depth))
        print('')
        print("SrF scattering rate: {:.3f} 1/s".format(SrF_scattering_rate))
        print("SrF heating rate: {:.3f} nK/s".format(SrF_heating_rate))
        print("SrF heating rate from 1064 nm laser: {:.3f} nK/s".format(SrF_shift_1064_final.heating_rate_nK_s))
        print("SrF heating rate from laser2 laser: {:.3f} nK/s".format(SrF_shift_laser2_final.heating_rate_nK_s))
        print('')
        print("Rb Stark shift: {:.2f} uk".format(Rb_shift_final))
        print("Rb trap depth w/ gravity: {:.2f} uk".format(Rb_trap_depth))
        print("Rb Stark shift from 1064 nm laser: {:.2f} uk".format(Rb_shift_1064_final.StarkShift_dict["scalar shift uK"]))
        print("Rb Stark shift from laser2: {:.2f} uk".format(Rb_shift_laser2_final.StarkShift_dict["scalar shift uK"]))
        print('')
        print("Rb scattering rate: {:.3f} 1/s".format(Rb_scattering_rate))
        print("Rb scattering rate from 1064 nm laser: {:.3f} 1/s".format(Rb_shift_1064_final.scattering_rate_invs))
        print("Rb scattering rate from laser2: {:.3f} 1/s".format(Rb_shift_laser2_final.scattering_rate_invs))
        print('')
        print("Rb heating rate: {:.3f} nK/s".format(Rb_heating_rate))
        print("Rb heating rate from 1064 nm laser: {:.3f} nK/s".format(Rb_shift_1064_final.heating_rate_nK_s))
        print("Rb heating rate from laser2 laser: {:.3f} nK/s".format(Rb_shift_laser2_final.heating_rate_nK_s))
        print('')
        print("1064 nm laser power: {:.4f} W".format(power_1064_final))
        print("laser2 ({:.2f} nm) power: {:.4f} W".format(wavelength_laser2, power_laser2_final))
        # print(Rb_shift_laser2_final.StarkShift_dict["scalar shift uK"])

    return {"SrF peak Stark shift ($\mu$K)":SrF_shift_final, "SrF trap depth w/ gravity ($\mu$K)":SrF_trap_depth, "SrF scattering rate (1/s)":SrF_scattering_rate, "SrF heating rate (nK/s)":SrF_heating_rate,
             "Rb peak Stark shift ($\mu$K)":Rb_shift_final, "Rb trap depth w/ gravity ($\mu$K)":Rb_trap_depth, "Rb scattering rate (1/s)":Rb_scattering_rate, "Rb heating rate (nK/s)":Rb_heating_rate,
            "1064 nm laser power (W)":power_1064_final, "laser2 power (W)":power_laser2_final, "laser2 wavelength (nm)":wavelength_laser2}

target_shift_Rb = 250 # uK
target_shift_SrF = 250 # uK
radius = 40 # 1/e^2 radius of laser beam, in um

headers = ["laser2 wavelength (nm)", "1064 nm laser power (W)", "laser2 power (W)", 
            "SrF peak Stark shift ($\mu$K)", "SrF trap depth w/ gravity ($\mu$K)", "SrF scattering rate (1/s)", "SrF heating rate (nK/s)", 
            "Rb peak Stark shift ($\mu$K)", "Rb trap depth w/ gravity ($\mu$K)", "Rb scattering rate (1/s)", "Rb heating rate (nK/s)"]
table = np.array([])
t0 = time.time()
for wavelength_laser2 in range(690, 761, 10):
    d = BichromaticTrap(target_shift_Rb, target_shift_SrF, wavelength_laser2, radius, print_result=False)
    for key in headers:
        val = "{:#.3g}".format(d[key])
        val = val[0:-1] if val[-1] == "." else val # remove decimal point at the end
        val = "\\textbf{"+f"{val}"+"}" if key in ["SrF peak Stark shift ($\mu$K)", "Rb peak Stark shift ($\mu$K)"] else val
        table = np.append(table, val)

    print("Wavelength {:d} nm calculation finished. (Elapsed time {:.1f} s.)".format(wavelength_laser2, time.time()-t0))

table = table.reshape((-1, len(headers)))
df = pd.DataFrame(table, columns=headers)
styler = df.style.hide(axis='index')
t = styler.to_latex(hrules=True, column_format=">{\centering}m{4.5em}"*(len(headers)-1)+">{\centering\\arraybackslash}m{4.5em}")

dir_path = os.path.dirname(os.path.realpath(__file__)) #+ "/BichromaticTrapTables/"
if not os.path.exists(dir_path):
    os.mkdir(dir_path) # make directory if it doesn't exist
with open(dir_path+f'BichromaticTrap_{target_shift_SrF}_{target_shift_Rb}.tex', 'w') as f:
    f.write(t)

print(t)