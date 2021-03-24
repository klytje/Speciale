import matplotlib.pyplot as plt
import sys
import math
from ROOT import TChain
from itertools import permutations
from dataframe import dataframe

chain = TChain('tree')
for f in input().split(): chain.AddFile(f)
fname = sys.argv[1]+"/"+sys.argv[0].split('/')[2].split('.')[0]

r = dataframe(chain)
r = r.define("esum","E_cm[0]+E_cm[1]+E_cm[2]")

xlim = (7.5e3,19e3)
ylim = (0,250e3)

plt.title(sys.argv[1])

plt.clf()
r.hist2d("exC12","p_tot",bins=400,range=(xlim,ylim),export=(fname,"count"),log=True)
plt.xlabel("$Ex({}^{12}\\mathrm{C})\\ [\\mathrm{keV}]$")
plt.ylabel("$p_\\mathrm{T}\\ [\\mathrm{keV}]$")
plt.gca().set_xlim(*xlim)
plt.gca().set_ylim(*ylim)

plt.savefig(fname+".pdf",dpi=600)

