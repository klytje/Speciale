import sys
import sympy
from sympy.physics.wigner import clebsch_gordan
from sympy.physics.wigner import racah
from sympy.functions.special.polynomials import legendre

# calculate b_n(l)
def calc_a(n, l):
	# we tell sympy to handle the division (S is short for sympify)
	return sympy.S(2*l*(l+1))/(2*l*(l+1) - n*(n+1))

# calculate the Biedenharn coefficient
def biedenharn(n, L, Lp, j, J):
	sign = sympy.Pow(-1, j-J-1)
	CG = clebsch_gordan(L, Lp, n, 1, -1, 0)
	W = racah(J, J, L, Lp, n, j)
	return sign*sympy.sqrt((2*L + 1)*(2*Lp + 1)*(2*J + 1)) * CG * W

# wrapper for non-mixed pure states
def calc_F(n, L, j, J):
	return biedenharn(n, L, L, j, J)

### main ###
# I've written this script with sympy, which symbolically solves equations. There are two reasons for this:
# 1: we get access to methods which calculate the clebsch-gordan and racah coefficients
# 2: we get symbolic results, which are nicer to work with (and can easily be converted to numerical values)
if len(sys.argv) != 3:
	print("Usage: python angular_correlation.py <j1> <l1> \nAssuming j = 2, j2 = 0, and l2 = 2")
	exit(1)

j1 = int(sys.argv[1]) # 12C
j = 2 # 8B
j2 = 0 # final alpha
l1 = int(sys.argv[2]) # first alpha (from 12C -> a + 8B)
l2 = 2 # second alpha (from 8B -> a + a)
print(f"Calculating the angular correlation for j1 = {j1}, j = {j}, j2 = {j2}, l1 = {l1}, l2 = {l2}\n")


W = sympy.S(0) # initialize W to 0
theta = sympy.symbols("\u03B8") # unicode for theta

# iterate over the summation symbol
for i in range(0, min([l1, l2, j])+1, 1):
	n = 2*i
	a1 = calc_a(n, l1)
	a2 = calc_a(n, l2)
	F1 = calc_F(n, l1, j1, j)
	F2 = calc_F(n, l2, j2, j)
	P = legendre(n, sympy.cos(theta))
	W += (a1*F1)*(a2*F2)*P

print(f"Symbolic result: W(\u03B8) = {W}")
print(f"Numeric result: W(\u03B8) = {W.evalf()}")