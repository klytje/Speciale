import sys
import mpmath
import numpy as np
import scipy.constants as const
from matplotlib import pyplot as plt

### units ###
hbar = const.hbar
c = const.c
fm = const.femto
MeV = 1e6*const.eV
m_p = 931.494*MeV/c**2

# hbar = 197.329
# c = 1
# fm = 1
# MeV = 1
# m_p = 931.494*MeV/c**2

### setup ###
# check input
if not (len(sys.argv) == 2):
	print("Usage: python width_converter <reduced width amplitude in keV^(1/2)>")
	exit(0)

# read input
gamma = float(sys.argv[1])*mpmath.sqrt(MeV/1000) # reduced width amplitude to be converted. input is in keV^(1/2), converted to MeV^(1/2)

# both particles are 4He
Z1 = Z2 = 2
A1 = A2 = 4

E_res = 3.03*MeV # resonance energy, 3030 keV
l = 2 # l quantum number of the excited state of 8Be
a = 1.2*fm*(mpmath.power(A1, 1.0/3) + mpmath.power(A2, 1.0/3)) # interaction radius
mu = A1*A2/(A1 + A2)*m_p # reduced mass

### calculate Gamma ###
def sommerfeld(rho):
	return const.alpha*Z1*Z2*c*mu*a/(rho*hbar)

# calculates the wavenumber based on the energy
def k(E):
	return mpmath.sqrt(2*mu*E)/hbar

# the regular coulomb wave function
def F(l, eta, rho):
	return mpmath.coulombf(l, eta, rho)

# the irregular coulomb wave function
def G(l, eta, rho):
	return mpmath.coulombg(l, eta, rho)

# the penetrability factor
def pen(l, rho):
	eta =  sommerfeld(rho)
	f = F(l, eta, rho)
	g = G(l, eta, rho)
	return rho/(f**2 + g**2)

# the shift factor
def shift(l, rho):
	eta =  sommerfeld(rho)
	f = F(l, eta, rho)
	g = G(l, eta, rho)

	# mpmath.diff requires that the function depends only on a single variable
	F_rho = lambda rho : F(l, eta, rho)
	G_rho = lambda rho : G(l, eta, rho)
	fp = mpmath.diff(F_rho, rho) # df/drho
	gp = mpmath.diff(G_rho, rho) # dg/drho

	return rho*(f*fp + g*gp)/(f**2 + g**2)

rho = a*k(E_res)
shift_rho = lambda rho : shift(l, rho)
Sp = mpmath.diff(shift_rho, rho)*rho/(2*E_res) # dS/dE = dS/drho * drho/dE = dS/drho * rho/2E
P = pen(l, rho)
S = shift(l, rho)

Gamma = 2*gamma**2*P/(1 + gamma**2*Sp*hbar)/MeV
print(f"The ordinary width corresponding to the reduced width {sys.argv[1]} keV^½ is")
print(f"{Gamma} MeV")

# Replicating the German stuff
x = np.linspace(0.01, 10, 50)*MeV # linspace of rho
# Z1 = 2
# Z2 = 4
# A1 = 4
# A2 = 8
# a = 5*fm
# mu = A1*A2/(A1 + A2)*m_p # reduced mass

# S0, S2 = [], []
# P0, P2 = [], []
# for i in x:
# 	rho = a*k(i)
# 	P0.append(float(pen(0, rho)))
# 	P2.append(float(pen(2, rho)))
# 	S0.append(float(shift(0, rho)))
# 	S2.append(float(shift(2, rho)))

# plt.figure()
# plt.title("Penetrability factor")
# plt.xlabel("E [MeV]")
# plt.ylabel("Pl")
# plt.plot(x/MeV, P0, label="l = 0")
# plt.plot(x/MeV, P2, label="l = 2")

# plt.figure()
# plt.title("Shift factor")
# plt.xlabel("E [MeV]")
# plt.ylabel("Sl")
# plt.plot(x/MeV, S0, label="l = 0")
# plt.plot(x/MeV, S2, label="l = 2")
# plt.show()


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
Z1 = 0 
Z2 = 0
x = np.linspace(0.01, 2, 50) # linspace of rho
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