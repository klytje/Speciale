import matplotlib.pyplot as plt
from ROOT import TChain
from dataframe import dataframe
from datatools import *
from mcfit import mcfitter
import numpy as np
from turbo import *
import matplotlib.gridspec as gridspec

np.seterr(all='ignore')

width =  1.0/72.27*393
plt.figure(figsize=(width,width/1.61803398875))


f = "interference"
datafolder = "data_tdc_cut"
labels = ["$E_\\mathrm{CM}\\ [\\mathrm{keV}]$","$\\varphi$","$\\rho$"]

parameters = ["eCM","phi","rho"]
ratios = [0,1]
legendof = {0: "$\\ell=1$", 1: "$\\ell=3$"}

plt.figure(figsize=(width,20*width/15))#/1.61803398875))

def plot_param(i):
    plt.subplot(3,1,i+1)
    for r in ratios:
        data,bins = load_data(datafolder)
        mcf = mcfitter(bins,data)
        balamuth,w = get_weighted(f)
        mcf.add_mc(balamuth,weight=w(r,0.5),label=f)
        mcf.fit(report=False)

        data,bins = load_data(datafolder,figure=parameters[i])
        mcf_ = mcfitter(bins,data)
        balamuth,w = get_weighted(f,figure=parameters[i])
        mcf_.add_mc(balamuth,weight=w(r,0.5),label=f)
        mcf_.val = np.array([1])
        #print(mcf_.chisquare())
        plt.hist(mcf_.bins[:-1],mcf_.bins,weights=mcf_.prediction(),histtype='step',label=legendof[r])
    plt.hist(mcf_.bins[:-1],mcf_.bins,weights=mcf_.data,histtype='step',label="Data",color='k')
    plt.xlabel(labels[i])
    plt.ylabel("Count")
    plt.xlim((min(bins),max(bins)))
    plt.legend(loc='upper left')

#for p in parameters:
    #plot_param(p,[0,1])
    #plt.legend()
    #plt.savefig("reportfig/model2/"+p+".pdf")
for i in range(3): plot_param(i)
plt.savefig("reportfig/model2/model2.pdf")
plt.savefig("reportfig/model2/model2.pgf")


def plot_d(chi2=False):
    plt.clf()
    #width =  1.0/72.27*393
    width =  1.0/72.27*413
    plt.figure(figsize=(width,width/1.61803398875*0.7))

    spec = gridspec.GridSpec(ncols=41, nrows=1)
    axc = plt.gcf().add_subplot(spec[0, 40])
    ax2 = plt.gcf().add_subplot(spec[0, 20:40])
    ax1 = plt.gcf().add_subplot(spec[0, 0:20],sharey=ax2)
    
    Hs = []
    for r in range(2):
        data,bins = load_data(datafolder)
        mcf = mcfitter(bins,data)
        balamuth,w = get_weighted(f)
        mcf.add_mc(balamuth,weight=w(r,0),label=f)
        mcf.fit()
        mcf.val = np.array([1])
        mcf.cov = mcf.cov
        print("dalitz_forfit",r,mcf.chisquare())

        data,bins = load_data(datafolder,figure="dalitzfig")
        mcf = mcfitter(bins,data)
        balamuth,w = get_weighted(f,figure="dalitzfig")
        mcf.add_mc(balamuth,weight=w(r,0),label=f)
        mcf.fit(report=False)
        mcf.val = np.array([1])
        mcf.cov = mcf.cov

        H = np.reshape(mcf.prediction(),(200,200))
        if chi2: H = np.reshape(-np.sign(mcf.data-mcf.prediction())*mcf.chisquare_binwise(),(200,200))
        Hs = Hs+[H]
    #H = np.reshape(mcf.data-mcf.prediction(),(200,200))

    mpl_data = RGBToPyCmap(turbo_colormap_data)
    plt.register_cmap(name='turbo', data=mpl_data, lut=turbo_colormap_data.shape[0])
    norm = [0,0]
    for H in Hs:
        norm[0] = min(np.min(H),norm[0])
        norm[1] = max(np.max(H),norm[1])
        if chi2: norm[1]=max(abs(norm[0]),abs(norm[1]))
    for H,ax in [(Hs[0],ax1),(Hs[1],ax2)]:
        if not chi2: H[H==0] = np.nan
        xedges = mcf.bins[0]
        yedges = mcf.bins[1]
        #print("r:",legendof[r],mcf.chisquare_binwise().sum())

        if not chi2: c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',cmap='turbo',vmin=norm[0],vmax=norm[1])
        else:c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',vmin=-350,vmax=350,cmap='seismic')
    #c = plt.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',vmin=-2e3,vmax=2e3,cmap='seismic')
        ax.set_aspect("equal")
        ax.axis('off')
    ax.set_ylim([-0.93,0.93])
    #plt.xlabel("$X$")
    #plt.ylabel("$Y$")
    #plt.colorbar(c)
    plt.colorbar(c,fraction=0.05,shrink=0.9,cax=axc)
    plt.sca(ax1)
    plt.title("$\\ell=1$")
    plt.sca(ax2)
    plt.title("$\\ell=3$")
    #lstr = "l1" if r==0 else "l3"
    chistr = "_chi2" if chi2 else ""
    plt.savefig("reportfig/model2/dalitz"+chistr+".pdf",dpi=1200)
    plt.savefig("reportfig/model2/dalitz"+chistr+".pgf",dpi=1200)

#plot_d(1)
plot_d()
#plot_d(1,True)
plot_d(True)
