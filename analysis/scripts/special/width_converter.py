import sys
import mpmath
import numpy as np
import scipy.constants as const
from matplotlib import pyplot as plt

# input
if not (len(sys.argv) == 2): 
    print("Usage: python width_converter <reduced width amplitude in keV^(1/2)>")

gamma = float(sys.argv[1]) * mpmath.sqrt(1000*const.eV) # reduced width amplitude to be converted

### units ###
hbar = 197.329 #hbar*c [MeV*fm]
hbar = 1
epsilon_0 = 1
e = mpmath.sqrt(4*const.pi*const.alpha)
m_p = 931.494 # measured in MeV
fm = 1

### setup parameters ###
# 4He
Z1 = 2
A1 = 2 

# 12C
Z2 = 2
A2 = 2

E_res = 3.03 # resonance energy, 3030 keV
l = 2 # l quantum number of the excited state of 8Be
a_c = 1.4*fm*(mpmath.power(A1, 1.0/3) + mpmath.power(A2, 1.0/3)) # interaction radius
M = A1*A2/(A1 + A2)*m_p # reduced mass

### calculate Gamma ###
def sommerfeld(rho):
    return Z1*Z2*M*a_c/(137.036*rho)

# momentum from beam energy
def k(E): 
    return mpmath.sqrt(2*M*E)/hbar

# the regular coulomb wave function
def F(l, eta, rho):
    return mpmath.coulombf(l, eta, rho)

# the irregular coulomb wave function
def G(l, eta, rho):
    return mpmath.coulombg(l, eta, rho)

# the penetrability factor
def pen(l, eta, rho): 
    f = F(l, eta, rho)
    g = G(l, eta, rho)
    return rho/(mpmath.absmax(f)*f + mpmath.absmax(g)*g)

# the shift factor
def shift(l, eta, rho):
    f = F(l, eta, rho)
    g = G(l, eta, rho)

    # mpmath.diff requires that the function depends only on a single variable
    F_rho = lambda rho : F(l, eta, rho) 
    G_rho = lambda rho : G(l, eta, rho) 
    fp = mpmath.diff(F_rho, rho) # df/drho
    gp = mpmath.diff(G_rho, rho) # dg/drho

    return pen(l, eta, rho)*(f*fp + g*gp)

rho = a_c*k(E_res) # E_res or E_beam ?
eta =  sommerfeld(rho) # the sommerfeld parameter
shift_rho = lambda rho : shift(l, eta, rho)
Sp = mpmath.diff(shift_rho, rho)*a_c*hbar*(-mpmath.power(E_res, -3.0/2)/(4*M)) # dS/dE
P = pen(l, eta, rho)

print(gamma**2*P)
print(gamma**2*Sp)

Gamma = 2*gamma**2*P/(1 + gamma**2*Sp)/1000 # convert to keV
print(f"The ordinary width corresponding to the reduced width {sys.argv[1]} keV^½ is")
print(f"{Gamma} keV")

x = np.linspace(0.01, 10, 50) # linspace of rho
# t1 = list(F(0, 0, z) for z in x)
# t2 = list(F(1, 0, z) for z in x)
# t3 = list(F(2, 0, z) for z in x)
# t4 = list(F(3, 0, z) for z in x)
# t5 = list(F(4, 0, z) for z in x)

# v1 = list(G(0, 0, z) for z in x)
# v2 = list(G(1, 0, z) for z in x)
# v3 = list(G(2, 0, z) for z in x)
# v4 = list(G(3, 0, z) for z in x)
# v5 = list(G(4, 0, z) for z in x)

# plt.figure()
# plt.title("Regular Coulomb functions")
# plt.xlabel("rho")
# plt.ylabel("Fl")
# plt.plot(x, t1, label="l = 0")
# plt.plot(x, t2, label="l = 1")
# plt.plot(x, t3, label="l = 2")
# plt.plot(x, t4, label="l = 3")
# plt.plot(x, t5, label="l = 4")
# plt.legend()

# plt.figure()
# plt.title("Irregular Coulomb functions")
# plt.xlabel("rho")
# plt.ylabel("Gl")
# plt.plot(x, v1, label="l = 0")
# plt.plot(x, v2, label="l = 1")
# plt.plot(x, v3, label="l = 2")
# plt.plot(x, v4, label="l = 3")
# plt.plot(x, v5, label="l = 4")
# plt.axis([0, 10, -2, 4])
# plt.legend()
# plt.show()

x = np.linspace(0.01, 2, 50) # linspace of rho
y1 = list(pen(0, 0, z) for z in x)
y2 = list(pen(1, 0, z) for z in x)
y3 = list(pen(2, 0, z) for z in x)
y4 = list(pen(3, 0, z) for z in x)
y5 = list(pen(4, 0, z) for z in x)

z1 = list(shift(0, 0, z) for z in x)
z2 = list(shift(1, 0, z) for z in x)
z3 = list(shift(2, 0, z) for z in x)
z4 = list(shift(3, 0, z) for z in x)
z5 = list(shift(4, 0, z) for z in x)

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