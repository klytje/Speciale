import os
import matplotlib.gridspec as gridspec
from turbo import *
import sys
from datatools import *
from anycache import anycache
import numpy as np
from mcfit import mcfitter
from iminuit import Minuit
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D

degree_folder = "/mnt/ssd/automake/"+sys.argv[1].split('/')[-1]
degree_folder = "reportfig/rmatrix"

data_folder = "data_tdc_cut"

width =  1.0/72.27*393
plt.figure(figsize=(width,width/1.61803398875))
plt.subplot(2,1,1)

folder = degree_folder+"/auto_8_1.42_36_2985_2-"

data,bins = load_data(data_folder)
balamuth,w = get_weighted("",path=folder)

def t(p,mc_weight):
    mc_weight[0] = w(p[0],p[1])
    return np.array([1]),mc_weight

mcf_ = mcfitter(bins,data,transform=t)
mcf_.add_mc(balamuth,weight=balamuth,label="Interference")
mcf_.fit(start=[0.2,0.7],limit=[(0,1)]*2,error=[0.4]*2)

plt.clf()
width =  1.0/72.27*413
plt.figure(figsize=(width,width/1.61803398875*0.7))
spec = gridspec.GridSpec(ncols=40, nrows=1)
ax2 = plt.gcf().add_subplot(spec[0, 21:40])
ax1 = plt.gcf().add_subplot(spec[0, 0:19],sharey=ax2)

def plot(ax,chi2=False):
    plt.sca(ax)
    data,bins = load_data(data_folder,figure="dalitzfig")
    balamuth,w = get_weighted("",figure="dalitzfig",path=folder)
    def t(p,mc_weight):
        mc_weight[0] = w(p[0],p[1])
        return np.array([1]),mc_weight
    mcf = mcfitter(bins,data,transform=t)
    mcf.add_mc(balamuth,weight=w(1,0))
    mcf.val = mcf_.val
    mcf.cov = mcf_.cov
    
    H = np.reshape(mcf.prediction(),(200,200))
    if chi2: H = np.reshape(-np.sign(mcf.data-mcf.prediction())*mcf.chisquare_binwise(),(200,200))

    mpl_data = RGBToPyCmap(turbo_colormap_data)
    plt.register_cmap(name='turbo', data=mpl_data, lut=turbo_colormap_data.shape[0])

    norm = [0,0]
    norm[0] = min(np.min(H),norm[0])
    norm[1] = max(np.max(H),norm[1])
    if chi2: norm[1]=max(abs(norm[0]),abs(norm[1]))
    if not chi2:H[H==0] = np.nan
    xedges = mcf.bins[0]
    yedges = mcf.bins[1]

    #c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',cmap='turbo',vmin=norm[0],vmax=norm[1])
    if not chi2: c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',cmap='turbo',vmin=norm[0],vmax=norm[1])
    else:c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',vmin=-200,vmax=200,cmap='seismic')
    plt.colorbar(c,fraction=0.05,shrink=1)

plot(ax1)
plot(ax2,chi2=True)
plt.ylim([-0.93,0.93])
for ax in (ax1,ax2):
    ax.axis('off')
    ax.set_aspect("equal")
plt.savefig("reportfig/rmatrix/rmatrix.pdf",dpi=600)
plt.savefig("reportfig/rmatrix/rmatrix.pgf",dpi=1200)
