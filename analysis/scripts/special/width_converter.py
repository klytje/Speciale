import sys
import sympy
import mpmath
import numpy as np
import scipy.constants as const
from matplotlib import pyplot as plt

# input
if not (len(sys.argv) == 2): 
    print("Usage: python width_converter <reduced width amplitude in keV^(1/2)>")

gamma = float(sys.argv[1]) * mpmath.sqrt(1000*const.eV) # reduced width amplitude to be converted

### setup parameters ###
# 4He
Z1 = 2
A1 = 2 
m1 = 4.0026*const.u

# 8Be
Z2 = 4
A2 = 4
m2 = 8.0053*const.u

E_res = 3030*1000*const.eV # resonance energy, 3030 keV
E_beam = 2001*1000*const.eV # beam energy, 2001 keV
l = 2 # l quantum number of the excited state of 8Be
a_c = 1.4*const.femto*(mpmath.power(A1, 1.0/3) + mpmath.power(A2, 1.0/3)) # interaction radius
M = m1*m2/(m1 + m2) # reduced mass

### calculate Gamma ###
# definition of the sommerfeld parameter
def sommerfeld(rho):
    return Z1*Z2*const.e**2*M*a_c/(4*const.pi*const.epsilon_0*const.hbar**2*rho)

# momentum from beam energy
def k(E): 
    return mpmath.sqrt(2*M*E)/const.hbar

# the regular coulomb wave function
def F(l, eta, rho):
    eta = sommerfeld(rho)
    return const.hbar**2/(2*M)*mpmath.coulombf(l, eta, rho)

# the irregular coulomb wave function
def G(l, eta, rho):
    return const.hbar**2/(2*M)*mpmath.coulombg(l, eta, rho)

# the penetrability factor
def pen(l, rho): 
    eta = sommerfeld(rho)
    f = F(l, eta, rho)
    g = G(l, eta, rho)
    return rho/(f**2 + g**2)

# the shift factor
def shift(l, rho):
    eta = sommerfeld(rho)
    f = F(l, eta, rho)
    g = G(l, eta, rho)

    # mpmath.diff requires that the function depends only on a single variable
    F_rho = lambda rho : F(l, eta, rho) 
    G_rho = lambda rho : G(l, eta, rho) 
    fp = mpmath.diff(F_rho, rho) # df/drho
    gp = mpmath.diff(G_rho, rho) # dg/drho

    return pen(l, rho)*(f*fp + g*gp)

rho = a_c*k(E_res) # E_res or E_beam ?
shift_rho = lambda rho : shift(l, rho)
Sp = mpmath.diff(shift_rho, rho)*a_c*const.hbar*(-mpmath.power(E_res, -3.0/2)/(4*M)) # dS/dE
P = pen(l, rho)

print(gamma**2*P)
print(gamma**2*Sp)

Gamma = 2*gamma**2*P/(1 + gamma**2*Sp) / (1000*const.eV) # convert to keV
print(f"The ordinary width corresponding to the reduced width {sys.argv[1]} keV^½ is")
print(f"{Gamma} keV")

x = np.linspace(0.01, 2, 50)
y1 = list(pen(0, z) for z in x)
y2 = list(pen(1, z) for z in x)
y3 = list(pen(2, z) for z in x)
y4 = list(pen(3, z) for z in x)
y5 = list(pen(4, z) for z in x)

z1 = list(shift(0, z) for z in x)
z2 = list(shift(1, z) for z in x)
z3 = list(shift(2, z) for z in x)
z4 = list(shift(3, z) for z in x)
z5 = list(shift(4, z) for z in x)

plt.figure()
plt.title("Penetration factor")
plt.xlabel("rho")
plt.ylabel("Pl")
plt.plot(x, y1, label="l = 0")
plt.plot(x, y2, label="l = 1")
plt.plot(x, y3, label="l = 2")
plt.plot(x, y4, label="l = 3")
plt.plot(x, y5, label="l = 4")
plt.legend()

plt.figure()
plt.title("Shift factors")
plt.xlabel("rho")
plt.ylabel("Sl")
plt.plot(x, z1, label="l = 0")
plt.plot(x, z2, label="l = 1")
plt.plot(x, z3, label="l = 2")
plt.plot(x, z4, label="l = 3")
plt.plot(x, z5, label="l = 4")
plt.legend()
plt.show()