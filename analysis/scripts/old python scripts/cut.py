import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
import math
from ROOT import TChain
from itertools import permutations
from dataframe import dataframe
import numpy as np
import matplotlib.gridspec as gridspec

width =  1.0/72.27*383
plt.figure(figsize=(width,width/1.61803398875))
#f, (ax, ax2) = plt.subplots(1, 2, sharey=True,figsize=(width,width/1.61803398875))

xlim = (0,11)
ylim = (0,200)

spec = gridspec.GridSpec(ncols=30, nrows=20)
axc = plt.gcf().add_subplot(spec[0:20, 29])
ax2 = plt.gcf().add_subplot(spec[16:20, 0:29])
ax1 = plt.gcf().add_subplot(spec[0:15, 0:29],sharex=ax2)


chain = TChain('a')
chain.Add("output/24*.root")
r_data = dataframe(chain).define("esum","eCM[0]+eCM[1]+eCM[2]").filter("mul==3")

d = r_data.asNumpy(["esum","pT"])

def plot_d(d):
    H, xedges, yedges = np.histogram2d(d["pT"]/1e3,d["esum"]/1e3,range=(ylim,xlim),bins=(400,500))
    H[H == 0] = np.nan
    c = plt.imshow(H[::-1,:], extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),vmin=1,aspect='auto',norm=LogNorm(vmin=1))

    Q = 9.305
    plt.plot([Q-.25, Q-.25],[0,500],'r')
    plt.plot([Q+.25, Q+.25],[0,500],'r')
    plt.plot([0, 10e6],[35,35],'r')
    plt.gca().set_ylim(*ylim)
    plt.ylabel("$p_\\mathrm{T}\\ [\\mathrm{MeV}]$")
    return c
c = plot_d(d)
#ax1.axes.xaxis.set_ticklabels([])
plt.setp(ax1.get_xticklabels(), visible=False)



#plt.colorbar(c)

chain = TChain('a')
chain.Add("output/reac-pB11-uniform.*.root")
r_uni = dataframe(chain).define("esum","eCM[0]+eCM[1]+eCM[2]").filter("mul==3")

plt.sca(ax2)
d = r_uni.asNumpy(["esum","pT"])
plot_d(d)

plt.sca(ax2)
plt.xlabel("$ \sum E_i^{\\mathrm{CM}}\\ [\\mathrm{MeV}]$")
plt.gca().set_xlim(*xlim)
plt.gca().set_ylim((0,50))

#cbar_ax = plt.gcf().add_axes([0.85, 0.15, 0.05, 0.7])
plt.gcf().colorbar(c, cax=axc)

plt.savefig("reportfig/cut.pdf",dpi=1200)
plt.savefig("reportfig/cut.pgf",dpi=1200)


