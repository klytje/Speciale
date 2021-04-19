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

# these units define MeV = fm = 1 to avoid precision errors
# hbar = 197.329
# c = 1
# fm = 1
# MeV = 1
# m_p = 931.494*MeV/c**2

### setup ###
# check input
if not (len(sys.argv) == 2 or len(sys.argv) == 3):
	print("Usage: python width_converter <reduced width amplitude in keV^(1/2)>")
	print("List mode outputs \"gamma [keV] Gamma [MeV]\" for easy piping (no units). List mode: -l")
	exit(0)

listmode = False
if (len(sys.argv) == 3): 
	if sys.argv[2] == "-l":
		listmode = True

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

Gamma = 2*gamma**2*P/(1 + gamma**2*Sp)/MeV
if listmode:
	print(f"{sys.argv[1]} {Gamma}")
else: 
	print(f"The ordinary width corresponding to the reduced width {sys.argv[1]} keV^½ is")
	print(f"{Gamma} MeV")

### debug ###
plot_german = False # creates a set of figures for 4He + 8Be which can be compared to the German thesis
plot_Baye = False # creates a plot of the shift and penetration factors for the neutral atom, which can be compared to the document by D. Baye

# replicates the German figures
if plot_german:
	E = np.linspace(0.01, 10, 50)*MeV # linspace of rho
	Z1 = 2
	Z2 = 4
	A1 = 4
	A2 = 8
	a = 5*fm
	mu = A1*A2/(A1 + A2)*m_p # reduced mass

	# penetrability and shift plots
	S0, S2 = [], []
	P0, P2 = [], []
	for e in E:
		rho = a*k(e)
		P0.append(float(pen(0, rho)))
		P2.append(float(pen(2, rho)))
		S0.append(float(shift(0, rho)))
		S2.append(float(shift(2, rho)))

	rho = a*k(E_res)
	B = shift(l, rho)
	dos = []
	# Be8 profile
	for e in E:
		rho = a*k(e)
		# shift_rho = lambda rho : shift(l, rho)
		# Sp = mpmath.diff(shift_rho, rho)*rho/(2*e) # dS/dE = dS/drho * drho/dE = dS/drho * rho/2E
		P = pen(l, rho)
		S = shift(l, rho)

		Gamma = 2*gamma**2*P
		Delta = -gamma**2*(S - B)

		dos.append(Gamma/((E_res + Delta - e)**2 + Gamma**2/4)) # density of states

	E /= MeV # set unit on x-axis
	# penetrability plot
	plt.figure()
	plt.title("Penetrability factor")
	plt.xlabel("E [MeV]")
	plt.ylabel(r"P_l")
	plt.plot(E, P0, label="l = 0")
	plt.plot(E, P2, label="l = 2")

	# shift plot
	plt.figure()
	plt.title("Shift factor")
	plt.xlabel("E [MeV]")
	plt.ylabel(r"S_l")
	plt.plot(E, S0, label="l = 0")
	plt.plot(E, S2, label="l = 2")

	# profile plot
	plt.figure()
	plt.title("8Be profile")
	plt.xlabel("E [MeV]")
	plt.ylabel("Density of states")
	plt.plot(E, dos)

	plt.show()

# replicates the neutral figure by Baye
if plot_Baye:
	Z1 = 0 
	Z2 = 0
	E = np.linspace(0.001, 2, 50)*hbar**2/(2*mu*a**2)
	pen0, pen1, pen2, pen3, pen4 = [], [], [], [], []
	shift0, shift1, shift2, shift3, shift4 = [], [], [], [], []
	for e in E:
		rho = a*k(e)
		pen0.append(pen(0, rho))
		pen1.append(pen(1, rho))
		pen2.append(pen(2, rho))
		pen3.append(pen(3, rho))
		pen4.append(pen(4, rho))

		shift0.append(shift(0, rho))
		shift1.append(shift(1, rho))
		shift2.append(shift(2, rho))
		shift3.append(shift(3, rho))
		shift4.append(shift(4, rho))

	E /= hbar**2/(2*mu*a**2)
	plt.figure()
	plt.title("Penetration factor")
	plt.xlabel(r"$E\ [\hbar^2 / 2\mu a^2]$")
	plt.ylabel(r"$P_l$")
	plt.plot(E, pen0, label="l = 0")
	plt.plot(E, pen1, label="l = 1")
	plt.plot(E, pen2, label="l = 2")
	plt.plot(E, pen3, label="l = 3")
	plt.plot(E, pen4, label="l = 4")
	plt.legend()

	plt.figure()
	plt.title("Shift factors")
	plt.xlabel(r"$E\ [\hbar^2 / 2\mu a^2]$")
	plt.ylabel(r"$S_l$")
	plt.plot(E, shift0, label="l = 0")
	plt.plot(E, shift1, label="l = 1")
	plt.plot(E, shift2, label="l = 2")
	plt.plot(E, shift3, label="l = 3")
	plt.plot(E, shift4, label="l = 4")
	plt.legend()
	plt.show()