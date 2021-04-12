import sys
import sympy
import mpmath
import numpy as np
import scipy.constants as const
from matplotlib import pyplot as plt

# input
if not (len(sys.argv) == 2): 
    print("Usage: python width_converter <reduced width amplitude in keV>")

gamma = float(sys.argv[1])*1000 # reduced width amplitude to be converted

### setup parameters ###
Z1 = 2 # 4He
A1 = 2 
Z2 = 4 # 8Be
A2 = 4 
m1 = 4.0026*const.u
m2 = 8.0053*const.u

beam = 2001*1000*const.eV # 2001keV
l = 2 # l quantum number of the excited state of 8Be
a_c = 1.4*const.femto*(mpmath.power(A1, 1.0/3) + mpmath.power(A2, 1.0/3)) # interaction radius

### calculate Gamma ###
# definition of the sommerfeld parameter
def sommerfeld(Z1, Z2, m1, m2, E):
    reduced_mass = lambda m1, m2 : m1*m2/(m1 + m2)
    return const.alpha*Z1*Z2*mpmath.sqrt(reduced_mass(m1, m2)*const.c**2/(2*E))

# momentum from beam energy
def k(E): 
    return mpmath.sqrt((m1**2 * const.c**2) - E**2/const.c**2) # E**2 = (mc**2)**2 + (pc)**2
#    return mpmath.sqrt(2*m1*E)/const.hbar

# the regular coulomb wave function
def F(E):
    eta = sommerfeld(Z1, Z2, m1, m2, E)
    return mpmath.coulombf(l, eta, a_c)

# the irregular coulomb wave function
def G(E):
    eta = sommerfeld(Z1, Z2, m1, m2, E)
    return mpmath.coulombg(l, eta, a_c)

# the penetrability factor
def pen(E): 
    rho = a_c*k(E)
    f = F(E)
    g = G(E)
    return rho/(f**2 + g**2)

# the shift factor
def shift(E):
    f = F(E)
    g = G(E)
    fp = mpmath.diff(F, E) # df/dE
    gp = mpmath.diff(G, E) # dg/dE
    return pen(E)*(f*fp + g*gp)

S = shift(beam)
Sp = mpmath.diff(shift, beam) # dS/dE
P = pen(beam)

Gamma = 2 * gamma**2 * P/(1 + gamma**2 * Sp) / (1000*const.eV)

# x = np.linspace(100, 200000, 100)
# y = list(pen(z) for z in x)
# plt.figure()
# plt.plot(x, y)
# plt.show()

print(f"The ordinary width corresponding to the reduced width {sys.argv[1]} keV^½ is")
print(f"{Gamma} keV")