import os
from turbo import *
import sys
import sys
from datatools import *
from anycache import anycache
import numpy as np
from mcfit import mcfitter
from iminuit import Minuit
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
    
mpl_data = RGBToPyCmap(turbo_colormap_data)
plt.register_cmap(name='turbo', data=mpl_data, lut=turbo_colormap_data.shape[0])

degree_folder = "/mnt/ssd/automake/"+sys.argv[1].split('/')[-1]
degree_folder = "reportfig/rmatrix"

data_folder = "figures/data_tdc_cut"
data0 = "../"


def fix(a,balamuth):
    a = a/balamuth
    a[np.isnan(a)] = 0
    return a

@anycache(cachedir='plotmake_cache')
def eval(coulomb_radius="12",r0Be="1.42",gBe="32.787",ExBe="3030",JC="2-",dfolder=data_folder,delta0=False):
    folder = degree_folder+"/figures/auto_"+coulomb_radius+"_"+r0Be+"_"+gBe+"_"+ExBe+"_"+JC

    data,bins = load_data("",path=dfolder)
    balamuth,w = get_weighted("",path=folder)

    
    if delta0:
        def t(p,mc_weight):
            mc_weight[0] = w(p[0],0)
            return np.array([1]),mc_weight
        mcf = mcfitter(bins,data,transform=t)
        mcf.add_mc(balamuth,weight=balamuth,label="Interference")
        mcf.fit(start=[0.5],limit=[(0,1)],error=[0.4])
        return mcf.chisquare()
    else:
        def t(p,mc_weight):
            mc_weight[0] = w(p[0],p[1])
            return np.array([1]),mc_weight
        mcf = mcfitter(bins,data,transform=t)
        mcf.add_mc(balamuth,weight=balamuth,label="Interference")
        mcf.fit(start=[0.2,0.7],limit=[(0,1)]*2,error=[0.9]*2)
        return mcf.chisquare()
    
    #return mcf.lnL_binwise().sum()

def parse(x): return x.split("_")[1:]
my_list = os.listdir(degree_folder+'/figures')
my_list = list(filter(lambda x: "auto" in x,my_list))
my_list = [parse(x) for x in my_list]

x = []
y = []
z1 = []
z2 = []

factor = 40

for a in my_list:
    if a[4] == "2-":
        x.append(float(a[2]))
        y.append(float(a[3])/factor)
        z1.append(eval(a[0],a[1],a[2],a[3],a[4],"../5_degree/figures/data_tdc_cut"))
        z2.append(eval(a[0],a[1],a[2],a[3],a[4],"../5_degree/figures/data_tdc_cut",delta0=True))

print(x,y,z1)
print(x,y,z2)

triang = mtri.Triangulation(x, y)

#fig = plt.figure(figsize=(width,width/1.61803398875))
width =  1.0/72.27*393
fig = plt.figure(figsize=(width,width*1.5))#/1.61803398875))
#fig = plt.figure()

def plotz(z,ax):
    xlim = [min(x+[32.787]), max(x+[32.787])]
    ylim = [min(y+[3030/factor]), max(y+[3030/factor])]
    N_ix = 20
    N_iy = 40

    xi, yi = np.meshgrid(np.linspace(xlim[0], xlim[1], N_ix), np.linspace(ylim[0],ylim[1], N_iy))
    interp_cubic_geom = mtri.CubicTriInterpolator(triang, z, kind='geom')
    zi_cubic_geom = interp_cubic_geom(xi, yi)

    vmin = np.min(zi_cubic_geom)
    vmax = np.max(zi_cubic_geom)
    vmax = vmin/2 + vmax/2*1.1


    c = ax.contourf(xi, yi, zi_cubic_geom,50,cmap='turbo_r',levels=np.linspace(vmin,vmax,50))
    ax.scatter(x,y, marker='.', c="black", alpha=0.4)
    ax.plot([32.787],[3037/factor], marker='x', color="k",markersize=3,ls='none')
    #plt.plot([v1[0],v2[0]],[v1[1]-1,v2[1]],marker='x',color='k',markersize=3,ls='none')
    
    print(vmin,int(vmin/100),int(vmax/100))
    r = range(int(vmin/100),int(vmax/100),int((vmax-vmin)/100/5))
    cbar = plt.colorbar(c,fraction=0.05,ticks=[100*r for r in r],pad=0.01)
    #cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
    #cbar.ax.set_yticklabels(['Low', 'high'])  # horizontal colorbar

#plt.title("$\chi^2$")
    ax.set_ylabel('$E_0$ $[\\mathrm{keV}]$')
    plt.xlim(xlim)
    plt.ylim(ylim)
    ticks = plt.gca().get_yticks()*factor
    plt.gca().set_yticklabels(ticks)

ax = fig.add_subplot(3,1,1)
plotz(z1,ax)
plt.title("varied $\\delta$ (startguess $\\delta=0.7$)")
plt.setp(ax.get_xticklabels(), visible=False)
ax = fig.add_subplot(3,1,2)
plotz(z2,ax)
plt.title("fixed $\\delta=0$")
ax.set_xlabel('$\\gamma_2$ $[\\sqrt{\\mathrm{keV}}]$')
ax = fig.add_subplot(3,1,3)
plotz([min(a,b) for a,b in zip(z1,z2)],ax)
plt.title("minimum of the two above")
ax.set_xlabel('$\\gamma_2$ $[\\sqrt{\\mathrm{keV}}]$')
plt.savefig("surf.pdf")
plt.savefig("surf.pgf")

