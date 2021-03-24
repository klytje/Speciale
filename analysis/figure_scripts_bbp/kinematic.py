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
ylim = (0,6000)

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
    .filter("Y<0.31")
    
ghist2d(r,["theta_lab[0]","theta_lab[1]","theta_lab[2]"],["E_lab[0]","E_lab[1]","E_lab[2]"],200,(xlim,ylim),fname)

#if "wU" in r.getColumnNames():
    #print("weighted?")
    #r.hist2d(["theta_lab[0]","theta_lab[1]","theta_lab[2]"],["eLab[0]","eLab[1]","eLab[2]"],bins=200,range=(xlim,ylim),weights=["wU*(0.22*f[0][0]+(1-0.22)*f[0][1]+2*sqrt(0.22*(1-0.22))*(f[0][2]*cos(2*3.1415*0.69)+f[0][3]*sin(2*3.1415*0.69)))"]*3)
#else:
    #r.hist2d(["theta_lab[0]","theta_lab[1]","theta_lab[2]"],["eLab[0]","eLab[1]","eLab[2]"],bins=200,range=(xlim,ylim))

plt.xlabel("$\\theta_{\\mathrm{Lab}}$")
plt.ylabel("$E_{\\mathrm{Lab}}$")
plt.gca().set_xlim(*xlim)
plt.gca().set_ylim(*ylim)

plt.savefig(fname+".pdf",dpi=1200)

