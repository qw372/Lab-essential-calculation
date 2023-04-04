import numpy as np
import matplotlib.pyplot as plt

def Snell_reflection(n1, theta1, n2, theta2):
    """
    Snell's law of reflection

    n1: index of reflection of incident medium
    theta1: incident angle
    n2: index of reflection of transmitted medium
    theta2: transmitted angle
    """

    rs = ((n1*np.cos(theta1)-n2*np.cos(theta2))/(n1*np.cos(theta1)+n2*np.cos(theta2)))**2 # reflection of s-polarized light
    rp = ((n1*np.cos(theta2)-n2*np.cos(theta1))/(n1*np.cos(theta2)+n2*np.cos(theta1)))**2 # reflection of s-polarized light

    return (rs+rp)/2

def Snell_transmission(n1, theta1, n2, theta2):
    """
    Snell's law of transmission

    n1: index of reflection of incident medium
    theta1: incident angle
    n2: index of reflection of transmitted medium
    theta2: transmitted angle
    """

    return 1 - Snell_reflection(n1, theta1, n2, theta2)

def fraction(l1, l2, l3, d1, d2, theta):
    """
    l1: distance between molecule cloud to the bottom of the lightpipe
    l2: length of light pipe
    l3: distance between PMT and top of lightpipe
    d1: diameter of lightpipe
    d2: diameter of PMT detection area
    theta: azimuthal angle of incident light

    See John Barry thesis Appendix J for more details
    """

    n = 1.55 # Quartz index of refraction

    if theta >= np.pi/2 or theta <= 0:
        return 0

    dd1 = l1*np.tan(theta)
    if dd1 >= d1/2:
        return 0

    theta_quartz = np.arcsin(np.sin(theta)/n)
    dd2 = l2*np.tan(theta_quartz)
    num_reflection = int(np.ceil((dd2+dd1-d1/2)/d1)) # number of reflections in quartz tube
    dd3 = d1*num_reflection+d1/2-dd1-dd2
    dd4 = l3*np.tan(theta)
    dd3 = dd3-dd4
    if dd3 < (d1-d2)/2 or dd3 > (d1+d2)/2:
        # if light leaks through the gap between lightpipe and PMT
        return 0

    t1 = Snell_transmission(1, theta, n, theta_quartz) # transmission of coupling from vacuum to quartz
    t2 = Snell_transmission(n, theta_quartz, 1, theta) # transmission of coupling from vacuum to quartz
    if n*np.sin(np.pi/2-theta_quartz) >= 1:
        r = 1
    else:
        theta2 = np.arcsin(n*np.sin(np.pi/2-theta_quartz))
        r = Snell_reflection(n, np.pi/2-theta_quartz, 1, theta2)
    
    return t1*t2*(r**num_reflection)


l1 = 15 # mm, distance between molecule cloud to the bottom of the lightpipe
l2 = 300 # mm, length of light pipe
l3 = 15 # mm, distance between PMT and top of lightpipe
d1 = 25*0.75 # mm, diameter of lightpipe
d2 = 22 # mm, diameter of PMT detection area

theta_list = np.linspace(0, np.pi/2, 5000)
dtheta = theta_list[1] - theta_list[0]
print("Fraction of light that is coupled into lightpipe can reach PMT: {:.3f}.".format(np.sum([fraction(l1, l2, l3, d1, d2, theta)*np.sin(theta)*dtheta for theta in theta_list])*(2*np.pi)/(2*np.pi*(1-np.cos(np.arctan(d1/2/l1))))))

fraction_list = np.zeros(len(theta_list))
for i, theta in enumerate(theta_list):
    fraction_list[i] = fraction(l1, l2, l3, d1, d2, theta)

# plt.plot(theta_list, fraction_list)
# plt.show()

