import matplotlib.pyplot as plt
from ROOT import TChain
from dataframe import dataframe
from datatools import *
from mcfit import mcfitter
import numpy as np
from turbo import *
import matplotlib.gridspec as gridspec

width =  1.0/72.27*403


f1 = "rgenl1"
f3 = "rgen"
datafolder = "true_events"

f_list = [[f1],[f3]]
parameters = ["eCM","phi","rho"]
labels = ["$E_\\mathrm{CM}\\ [\\mathrm{keV}]$","$\\varphi$","$\\rho$"]

legendof = {f1: "$\\ell=1$", f3: "$\\ell=3$"}

plt.figure(figsize=(width,20*width/15))#/1.61803398875))

def plot_param(i):
    plt.subplot(3,1,i+1)
    for files in f_list:
        data,bins = load_data(datafolder,figure=parameters[i])
        mcf_ = mcfitter(bins,data)
        for f in files:
            m,bins = load_data(f,figure=parameters[i])
            mcf_.add_mc(m,label=f)
        mcf_.val = np.array([1])
        print(parameters[i],legendof[files[0]],mcf_.chisquare())

        plt.hist(mcf_.bins[:-1],mcf_.bins,weights=mcf_.prediction(),histtype='step',label=legendof[files[0]])
    plt.hist(mcf_.bins[:-1],mcf_.bins,weights=mcf_.data,histtype='step',label="Data",color='k')
    plt.xlabel(labels[i])
    plt.ylabel("Count")
    plt.xlim((min(bins),max(bins)))
    plt.legend(loc='upper left')

for i in range(3): plot_param(i)
plt.savefig("reportfig/model1/model1.pdf")
plt.savefig("reportfig/model1/model1.pgf")

def plot_d(chi2=False):
    plt.clf()
    width =  1.0/72.27*413
    plt.figure(figsize=(width,width/1.61803398875*0.7))

    spec = gridspec.GridSpec(ncols=41, nrows=1)
    axc = plt.gcf().add_subplot(spec[0, 40])
    ax2 = plt.gcf().add_subplot(spec[0, 20:40])
    ax1 = plt.gcf().add_subplot(spec[0, 0:20],sharey=ax2)

    Hs = []
    for f in [f1,f3]:
        data,bins = load_data(datafolder)
        mcf = mcfitter(bins,data)
        m,bins = load_data(f)
        mcf.add_mc(m,label=f)
        mcf.val = np.array([1])
        print("dalitz_forfit",f,mcf.chisquare())

        data,bins = load_data(datafolder,figure="dalitzfig")
        mcf = mcfitter(bins,data)
        m,bins = load_data(f,figure="dalitzfig")
        mcf.add_mc(m,label=f)
        mcf.val = np.array([1])
        print("dalitz",f,mcf.chisquare())
        H = np.reshape(mcf.prediction(),(200,200))
        if chi2: H = np.reshape(-np.sign(mcf.data-mcf.prediction())*mcf.chisquare_binwise(),(200,200))
        Hs = Hs+[H]

    norm = [0,0]
    for H in Hs:
        norm[0] = min(np.min(H),norm[0])
        norm[1] = max(np.max(H),norm[1])
        if chi2: norm[1]=max(abs(norm[0]),abs(norm[1]))
    
    mpl_data = RGBToPyCmap(turbo_colormap_data)
    plt.register_cmap(name='turbo', data=mpl_data, lut=turbo_colormap_data.shape[0])

    for H,ax in [(Hs[0],ax1),(Hs[1],ax2)]:
        H[H==0] = np.nan
        xedges = mcf.bins[0]
        yedges = mcf.bins[1]

        if not chi2: c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',cmap='turbo',vmin=norm[0],vmax=norm[1])
        else:c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',vmin=-400,vmax=400,cmap='seismic')
        ax.set_aspect("equal")
        ax.axis('off')
    ax.set_ylim([-0.93,0.93])
    plt.colorbar(c,fraction=0.05,shrink=0.9,cax=axc)
    
    plt.sca(ax1)
    plt.title("$\\ell=1$")
    plt.sca(ax2)
    plt.title("$\\ell=3$")
    
    plt.savefig("reportfig/model1/rgen"+("_chi2"if chi2 else "")+".pdf",dpi=1200)
    plt.savefig("reportfig/model1/rgen"+("_chi2"if chi2 else "")+".pgf",dpi=1200)

plot_d()
plot_d(True)
