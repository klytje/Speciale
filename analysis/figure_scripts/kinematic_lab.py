import matplotlib.pyplot as plt
import sys
import math
from ROOT import TChain
from itertools import permutations
from dataframe import *

chain = TChain('tree')
for f in input().split(): chain.AddFile(f)
fname = sys.argv[1]+"/"+sys.argv[0].split('/')[2].split('.')[0]

r = dataframe(chain)
r = r.define("esum","E_cm[0]+E_cm[1]+E_cm[2]")

xlim = (25,170)
ylim = (-180,180)

plt.figure(figsize=(6.4,6.4/1.85))
plt.title(sys.argv[1])
r = r.filter("p_tot<50e3")\
    .filter("abs(deltaE)<200")\
    .define("E2","min(E_cm[0],min(E_cm[1],E_cm[2]))")\
    .define("E1","max(min(E_cm[0],E_cm[1]), min(max(E_cm[0],E_cm[1]),E_cm[2]))")\
    .define("E0","max(E_cm[0],max(E_cm[1],E_cm[2]))")\
    .define("X","((E0+2*E1)/esum-1)/sqrt(3)")\
    .define("Y","E0/esum-1.0/3.0")\
    .filter("pow(X,2)+pow(Y,2)<1.0/9.0")\
    .filter("E_cm[0]<6000 && E_cm[1]<6000 && E_cm[2]<6000")
#r.hist2d(["theta_lab[0]","theta_lab[1]","theta_lab[2]"],["phiLab[0]","phiLab[1]","phiLab[2]"],bins=500,range=(xlim,ylim),export=(fname,"count"))
    #.filter("E_cm[0]>1000 && E_cm[1]>1000 && E_cm[2]>1000")\

ghist2d(r,["theta_lab[0]","theta_lab[1]","theta_lab[2]"],["phi_lab[0]","phi_lab[1]","phi_lab[2]"],500,(xlim,ylim),fname)

#plt.clf()
#if "wU" in r.getColumnNames():
    #r.hist2d(["theta_lab[0]","theta_lab[1]","theta_lab[2]"],["phiLab[0]","phiLab[1]","phiLab[2]"],bins=500,range=(xlim,ylim),weights=["wU*(0.22*f[0][0]+(1-0.22)*f[0][1]+2*sqrt(0.22*(1-0.22))*(f[0][2]*cos(2*3.1415*0.69)+f[0][3]*sin(2*3.1415*0.69)))"]*3)
#else:
    #r.hist2d(["theta_lab[0]","theta_lab[1]","theta_lab[2]"],["phiLab[0]","phiLab[1]","phiLab[2]"],bins=500,range=(xlim,ylim))

plt.xlabel("$\\theta_{\\mathrm{lab}}$")
plt.ylabel("$\\phi_{\\mathrm{lab}}$")
plt.gca().set_xlim(*xlim)
plt.gca().set_ylim(*ylim)
#plt.gca().set_aspect("equal")

plt.savefig(fname+".pdf",dpi=600)

