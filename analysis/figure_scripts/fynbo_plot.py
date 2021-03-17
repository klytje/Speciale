import matplotlib.pyplot as plt
import sys
import math
from ROOT import TChain
from itertools import permutations
from dataframe import dataframe

chain = TChain('a')
for f in input().split(): chain.AddFile(f)
fname = sys.argv[1]+"/"+sys.argv[0].split('/')[2].split('.')[0]

r = dataframe(chain)
r = r.define("esum","eCM[0]+eCM[1]+eCM[2]")

xlim = (0,12000)
ylim = (0,6500)

plt.title(sys.argv[1])
r = r.filter("pT<50e3")

plt.clf()
r.hist2d(["esum","esum","esum"],["eCM[0]","eCM[1]","eCM[2]"],bins=400,range=(xlim,ylim),export=(fname,"count"),log=True)

plt.xlabel("$E_{\\mathrm{tot}}\\ [\\mathrm{keV}]$")
plt.ylabel("$E_{\\alpha}\\ [\\mathrm{keV}]$")
plt.gca().set_xlim(*xlim)
plt.gca().set_ylim(*ylim)

plt.savefig(fname+".pdf",dpi=600)

