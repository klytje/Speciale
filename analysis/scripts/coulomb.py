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

degree_folder = "../autogen/"+sys.argv[1].split('/')[-1]
degree_folder = "reportfig/coulomb"

data_folder = "data_tdc_cut"

@anycache(cachedir="reportfig/coulomb/cache")
def eval(coulomb_radius="12",r0Be="1.42",gBe="32.787",ExBe="3030",JC="2-",figure="dalitz_forfit"):
    folder = degree_folder+"/figures/auto_"+coulomb_radius+"_"+r0Be+"_"+gBe+"_"+ExBe+"_"+JC

    data,bins = load_data(data_folder,figure=figure)
    balamuth,w = get_weighted("",path=folder,figure=figure)

    def t(p,mc_weight):
        mc_weight[0] = w(p[0],p[1])
        return np.array([1]),mc_weight

    mcf = mcfitter(bins,data,transform=t)
    mcf.add_mc(balamuth,weight=balamuth,label="Interference")
    mcf.fit(start=[0.2,0.7],limit=[(0,1)]*2,error=[0.4]*2)
    
    #return mcf.lnL_binwise().sum()
    return mcf.chisquare()

def parse(x): return x.split("_")[1:]
my_list = os.listdir(degree_folder+'/figures')
my_list = list(filter(lambda x: "auto" in x,my_list))
my_list = [parse(x) for x in my_list]

x = []
y = []
y2 = []

for a in my_list:
    x.append(float(a[0]))
    print(a[0])
    y.append(eval(a[0],a[1],a[2],a[3],a[4]))
    y2.append(eval(a[0],a[1],a[2],a[3],a[4],figure="eCM"))

print(x,y)

width =  1.0/72.27*393
plt.figure(figsize=(width,width/1.61803398875))
plt.subplot(2,1,1)

x_,y = zip(*sorted(zip(x,y)))
plt.plot(x_,y,label="Dalitz fit")
plt.plot([8,8],[0,1e6],'k',linestyle='--')
plt.xlim([min(x_),max(x_)])
plt.ylim(min(y)*0.98,y[0]/2+y[-1]/2)
plt.setp(plt.gca().get_xticklabels(), visible=False)
x_,y = zip(*sorted(zip(x,y2)))
plt.ylabel("$\\chi^2$")
plt.legend()

plt.subplot(2,1,2)
plt.plot([0,0],[0,0])
plt.plot(x_,y,label="$E_{\\mathrm{CM}}$ fit")
plt.plot([10,10],[0,1e6],'k',linestyle='--')
plt.xlim([min(x_),max(x_)])
plt.ylim(min(y)*0.98,y[0]/2+y[-1]/2)
plt.xlabel("Channel radius $a_c'$ $[\\mathrm{fm}]$")
plt.ylabel("$\\chi^2$")
plt.legend()

plt.savefig("reportfig/coulomb/coulomb.pdf")
plt.savefig("reportfig/coulomb/coulomb.pgf")




folder = degree_folder+"/auto_8_1.42_32.79_3037_2-"

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
    #H[7,:] = 0
    if not chi2:H[H==0] = np.nan
    xedges = mcf.bins[0]
    yedges = mcf.bins[1]

    #c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',cmap='turbo',vmin=norm[0],vmax=norm[1])
    if not chi2: c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',cmap='turbo',vmin=norm[0],vmax=norm[1])
    else:c = ax.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',vmin=-100,vmax=100,cmap='seismic')
    plt.colorbar(c,fraction=0.05,shrink=1)

plot(ax1)
plot(ax2,chi2=True)
plt.ylim([-0.93,0.93])
for ax in (ax1,ax2):
    ax.axis('off')
    ax.set_aspect("equal")
plt.savefig("reportfig/coulomb/coulomb_compare.pdf",dpi=600)
plt.savefig("reportfig/coulomb/coulomb_compare.pgf",dpi=1200)
