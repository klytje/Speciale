import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from ROOT import TChain
from dataframe import dataframe
from datatools import *
from mcfit import mcfitter
import numpy as np
from turbo import *

width =  1.0/72.27*393


files = {
        "delta0_coulomb_restricted": [42.5, 3000],
        "deltafree_coulomb_restricted": [36.5, 3020],
        "delta0_coulomb": [41.2, 2955],
        "delta0_nocoulomb": [40.1, 2800],
        "deltafree_coulomb": [36, 2975],
        "deltafree_nocoulomb": [37, 2850]
}


datafolder = "data_tdc_cut"


def plot_d(f):
    parsed = f.split(".")[0].split("_")
    deltafree = parsed[0]=="deltafree"
    coulomb = parsed[1]=="coulomb"
    restricted = len(parsed)==3
    fname = "_restringeret" if restricted else ""



    plt.clf()
    width =  1.0/72.27*413
    plt.figure(figsize=(width,width/1.61803398875*0.7))

    spec = gridspec.GridSpec(ncols=50, nrows=50)
    ax1 = plt.gcf().add_subplot(spec[0:50, 0:20])
    axc1 = plt.gcf().add_subplot(spec[3:47, 20:21])
    ax2 = plt.gcf().add_subplot(spec[0:500, 25:45],sharey=ax1)
    axc2 = plt.gcf().add_subplot(spec[3:47, 45:46])


    k = 0
    delta = 0
    chisquare = 0
    ndf = 0
    
    mcf = None
    if deltafree:
        starts=[0.23,0.67]
        def t(p,mc_weight):
            mc_weight[0] = w(p[0],p[1])
            return np.array([1]),mc_weight
        data,bins = load_data(datafolder,figure="dalitz_forfit"+fname)
        mcf_ = mcfitter(bins,data,transform=t)
        balamuth,w = get_weighted(f,figure="dalitz_forfit"+fname)
        mcf_.add_mc(balamuth,weight=balamuth,label=f)
        mcf_.fit(start=starts,error=[0.1]*len(starts),limit=[(0,1)]*len(starts))
        chisquare = mcf_.chisquare()
        ndf = mcf_.ndf()

        data,bins = load_data(datafolder,figure="dalitzfig"+fname)
        balamuth,w = get_weighted(f,figure="dalitzfig"+fname)
        mcf = mcfitter(bins,data,transform=t)
        mcf.add_mc(balamuth,weight=w(1,0),label=f)
        mcf.val = mcf_.val
        mcf.cov = mcf_.cov
        k = round(mcf.val[0],3)
        delta = round(mcf.val[1],3)
    else:
        starts=[0.5]
        def t(p,mc_weight):
            mc_weight[0] = w(p[0],0)
            return np.array([1]),mc_weight
        data,bins = load_data(datafolder,figure="dalitz_forfit"+fname)
        mcf_ = mcfitter(bins,data,transform=t)
        balamuth,w = get_weighted(f,figure="dalitz_forfit"+fname)
        mcf_.add_mc(balamuth,weight=balamuth,label=f)
        mcf_.fit(start=starts,error=[0.1]*len(starts),limit=[(0,1)]*len(starts))
        chisquare = mcf_.chisquare()
        ndf = mcf_.ndf()

        data,bins = load_data(datafolder,figure="dalitzfig"+fname)
        balamuth,w = get_weighted(f,figure="dalitzfig"+fname)
        mcf = mcfitter(bins,data,transform=t)
        mcf.add_mc(balamuth,weight=w(1,0),label=f)
        mcf.val = mcf_.val
        mcf.cov = mcf_.cov
        k = round(mcf.val[0],3)

    def plot(ax,axc,chi2=False,restricted=False):
        H = np.reshape(mcf.prediction(),(200,200))
        if chi2: H = np.reshape(-np.sign(mcf.data-mcf.prediction())*mcf.chisquare_binwise(),(200,200))

        norm = [0,0]
        norm[0] = min(np.min(H),norm[0])
        norm[1] = max(np.max(H),norm[1])
        if chi2: norm[1]=max(abs(norm[0]),abs(norm[1]))
        if not chi2:H[H==0] = np.nan
        xedges = mcf.bins[0]
        yedges = mcf.bins[1]

        mpl_data = RGBToPyCmap(turbo_colormap_data)
        plt.register_cmap(name='turbo', data=mpl_data, lut=turbo_colormap_data.shape[0])

        if not chi2: c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',cmap='turbo',vmin=norm[0],vmax=norm[1])
        else: c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',vmin=-norm[1],vmax=norm[1],cmap='seismic')
        
        #if restricted:
            #phi = np.arange(0,2*3.14159,0.1)
            #x = 0.55 + np.cos(phi)*0.27
            #y = 1 + np.sin(phi)*0.27
            #ax.plot(x,y,'-k')
            #ax.text(0.35,0.9,"removed",size="xx-small")
       
        plt.colorbar(c,fraction=0.01,shrink=0.1,cax=axc)
        ax.set_aspect("equal")
        ax.axis('off')
        ax.set_ylim([-0.93,0.93])



    
    plot(ax1,axc1,restricted=restricted)
    plot(ax2,axc2,chi2=True,restricted=restricted)
    coulombstring = "8\\mathrm{fm}" if coulomb else "\\mathrm{None}"
    title = "$\\gamma_2="+str(files[f][0])+"\\mathrm{keV}^{1/2}$, $E_0="+str(files[f][1])+"\\mathrm{keV}$, $a_c\'="+ coulombstring +"$, $k="+str(k)+"$,$\\delta="+str(delta)+"$, $\\chi^2/\\mathrm{ndf}="+str(round(chisquare))+"/"+str(round(ndf))+"$"
    plt.suptitle(title,size='xx-small')
    plt.savefig("articlefig/"+f+".pdf",dpi=1200)
    #plt.savefig("articlefig/dalitz"+chistr+".pgf",dpi=1200)

for f in files:
    plot_d(f)
