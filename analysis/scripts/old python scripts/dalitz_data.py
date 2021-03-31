import matplotlib.pyplot as plt
from ROOT import TChain
from dataframe import dataframe
from datatools import *
from mcfit import mcfitter
import numpy as np
import math
from turbo import *
import matplotlib.gridspec as gridspec

datafolder = "data_tdc_cut"
width =  1.0/72.27*413
plt.figure(figsize=(width/2,width/1.61803398875*0.7))

data,bins = load_data(datafolder,figure="dalitzfig")
mcf = mcfitter(bins,data)

H = np.reshape(mcf.data,(200,200))

mpl_data = RGBToPyCmap(turbo_colormap_data)
plt.register_cmap(name='turbo', data=mpl_data, lut=turbo_colormap_data.shape[0])
H[H==0] = np.nan
xedges = mcf.bins[0]
yedges = mcf.bins[1]

c = plt.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',cmap='turbo')
plt.gca().set_aspect("equal")
plt.gca().axis('off')
#scale = (1.0-250/(9325/3))/math.cos(math.pi/6)
# cos(d) X = l = 1-250keV/9300
#plt.plot([0,scale],[scale*math.tan(math.pi/3),0],'r')
plt.ylim([-1,0.93])
plt.colorbar(c,fraction=0.05,shrink=0.9)
plt.savefig("reportfig/dalitz_data.pdf",dpi=1200)
plt.savefig("reportfig/dalitz_data.pgf",dpi=1200)

