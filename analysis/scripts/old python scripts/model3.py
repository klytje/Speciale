import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from ROOT import TChain
from dataframe import dataframe
from datatools import *
from mcfit import mcfitter
import numpy as np
from turbo import *

width =  1.0/72.27*393
#plt.figure(figsize=(width,width/1.61803398875))

labels = ["$E_\\mathrm{CM}\\ [\\mathrm{keV}]$","$\\varphi$","$\\rho$"]

f = "interference"
datafolder = "data_tdc_cut"

parameters = ["eCM","phi","rho"]
starts=[[0.23,0.67],[0.6,0.05]]
starts = [[0.2,0.67],[0.1]]
#legendof = {starts[0]: "$\\ell=1$", starts[1]: "$\\ell=3$"}

plt.figure(figsize=(width,20*width/15))#/1.61803398875))

def plot_param(i):
    plt.subplot(3,1,i+1)
    for j in range(2):
        data,bins = load_data(datafolder)
        balamuth,w = get_weighted(f)
        def t(p,mc_weight):
            mc_weight[0] = w(p[0],p[1])
            return np.array([1]),mc_weight
        def t2(p,mc_weight):
            mc_weight[0] = p[0]*w(0.0,0.0) + (1-p[0])*w(1.0,0.0)
            return np.array([1]),mc_weight
        if j==0:
            mcf = mcfitter(bins,data,transform=t)
        else:
            mcf = mcfitter(bins,data,transform=t2)
        mcf.add_mc(balamuth,weight=balamuth,label=f)
        mcf.fit(start=starts[j],error=[0.5]*len(starts[j]),limit=[(0,1)]*len(starts[j]))

        data,bins = load_data(datafolder,figure=parameters[i])
        balamuth,w = get_weighted(f,figure=parameters[i])
        #def t(p,mc_weight):
            #mc_weight[0] = w(p[0],p[1])
            #return np.array([1]),mc_weight
        if j==0:
            mcf_ = mcfitter(bins,data,transform=t)
        else:
            mcf_ = mcfitter(bins,data,transform=t2)
        mcf_.add_mc(balamuth,weight=balamuth,label=f)
        mcf_.val = mcf.val
        mcf_.cov = mcf.cov
        print(mcf_.prediction().sum(),mcf_.data.sum())
        if j==0:
            label = "$k="+("%.2f"%(mcf.val[0]*100))+"\\%$\n$\\delta=2\\pi\\cdot"+("%.2f"%(mcf.val[1]*100))+"\\%$"
        else:
            label = "$k="+("%.2f"%(100-mcf.val[0]*100))+"\\%$"
        plt.hist(mcf_.bins[:-1],mcf_.bins,weights=mcf_.prediction(),histtype='step',label=label)
    plt.hist(mcf_.bins[:-1],mcf_.bins,weights=mcf_.data,histtype='step',label="Data",color='k')
    plt.xlabel(labels[i])
    plt.ylabel("Count")
    plt.xlim((min(bins),max(bins)))
    plt.legend(loc='upper left')

#for p in parameters:
    #plt.figure(figsize=(width,width/1.61803398875))
    #plot_param(p)
    #plt.legend()
    #plt.savefig("reportfig/model3/"+p+".pdf")
for i in range(3): plot_param(i)
plt.savefig("reportfig/model3/model3.pdf")
plt.savefig("reportfig/model3/model3.pgf")


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
    for j in range(2):
        if j==0:
            def t(p,mc_weight):
                mc_weight[0] = w(p[0],p[1])
                return np.array([1]),mc_weight
            data,bins = load_data(datafolder)
            mcf_ = mcfitter(bins,data,transform=t)
        else:
            def t(p,mc_weight):
                mc_weight[0] = p[0]*w(0.0,0.0) + (1-p[0])*w(1.0,0.0)
                return np.array([1]),mc_weight
            data,bins = load_data(datafolder)
            mcf_ = mcfitter(bins,data,transform=t)
        balamuth,w = get_weighted(f)
        mcf_.add_mc(balamuth,weight=balamuth,label=f)
        mcf_.fit(start=starts[j],error=[0.1]*len(starts[j]),limit=[(0,1)]*len(starts[j]))
    
        data,bins = load_data(datafolder,figure="dalitzfig")
        balamuth,w = get_weighted(f,figure="dalitzfig")
        mcf = mcfitter(bins,data,transform=t)
        mcf.add_mc(balamuth,weight=w(1,0),label=f)
        mcf.val = mcf_.val
        mcf.cov = mcf_.cov


        H = np.reshape(mcf.prediction(),(200,200))
        if chi2: H = np.reshape(-np.sign(mcf.data-mcf.prediction())*mcf.chisquare_binwise(),(200,200))
        Hs = Hs+[H]

    mpl_data = RGBToPyCmap(turbo_colormap_data)
    plt.register_cmap(name='turbo', data=mpl_data, lut=turbo_colormap_data.shape[0])
    norm = [0,0]
    for H in Hs:
        norm[0] = min(np.min(H),norm[0])
        norm[1] = max(np.max(H),norm[1])
        if chi2: norm[1]=max(abs(norm[0]),abs(norm[1]))
    for H,ax in [(Hs[0],ax1),(Hs[1],ax2)]:
        if not chi2:H[H==0] = np.nan
        xedges = mcf.bins[0]
        yedges = mcf.bins[1]

        if not chi2: c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',cmap='turbo',vmin=norm[0],vmax=norm[1])
        else:c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',vmin=-200,vmax=200,cmap='seismic')
    #c = plt.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',vmin=-2e3,vmax=2e3,cmap='seismic')
        ax.set_aspect("equal")
        ax.axis('off')
    ax.set_ylim([-0.93,0.93])
    #plt.xlabel("$X$")
    #plt.ylabel("$Y$")
    #plt.colorbar(c)
    plt.colorbar(c,fraction=0.05,shrink=0.9,cax=axc)
    plt.sca(ax1)
    plt.title("Interference")
    plt.sca(ax2)
    plt.title("No interference")
    #lstr = "l1" if r==0 else "l3"
    chistr = "_chi2" if chi2 else ""
    plt.savefig("reportfig/model3/dalitz"+chistr+".pdf",dpi=1200)
    plt.savefig("reportfig/model3/dalitz"+chistr+".pgf",dpi=1200)

plot_d()
plot_d(True)
