from figure_scripts.dataframe import dataframe
from itertools import permutations
from anycache import anycache
from ROOT import TChain
import matplotlib.pyplot as plt
import numpy as np
from iminuit import Minuit
import math

chain = TChain('m')
for f in input().split(): chain.AddFile(f)

r = dataframe(chain)
b = r

## eff
c = r
c.hist(["eLab[0]","eLab[1]","eLab[2]"],bins=200,range=(200,800),export=("tdc/eff/testname.npz","count"))
for i in range(3):
    c2 = c.filter("BT["+str(i)+"]==0")
    c2.hist("eLab["+str(i)+"]",bins=200,range=(200,800),export=("tdc/eff/testname"+str(i)+".npz","count"))

## eff end

FT = "FT_"
BT = "BT_"

for i in range(3):
    j = str(i)
    r = r.define(FT+j,"(FT["+j+"]<6.53e7) ? FT["+j+"]+409500 : FT["+j+"]")
    r = r.define(BT+j,"(BT["+j+"]<6.53e7) ? BT["+j+"]+409500 : BT["+j+"]")

r_back = r

#
# fig tdc_raw
#
for i in range(3):
    j = str(i)
    b = b.filter("FT["+j+"]>1000")
    b = b.filter("BT["+j+"]>1000")
    b = b.define(FT+j,"(FT["+j+"]%651264)/100")
    b = b.define(BT+j,"(BT["+j+"]%651264)/100")
plt.figure()
width =  1.0/72.27*393
plt.figure(figsize=(width,width/1.61803398875))
f, (ax, ax2) = plt.subplots(1, 2, sharey=True,figsize=(width,width/2))#1.61803398875))
lim = (0,4096)
plt.sca(ax)
ax.set_xlim((0,718))
ax.spines['right'].set_visible(False)
ax.yaxis.set_ticks_position('left')
b.hist2d(FT+"0",BT+"0",bins=(400,400),range=(lim,lim),log=True)
ax.images[-1].colorbar.remove()
plt.ylabel("TDC channel (back)")

plt.sca(ax2)
ax2.spines['left'].set_visible(False)
ax2.yaxis.set_ticks_position('right')
ax2.tick_params(labelright='off')
ax2.set_xlim((3521,4096))
ax2.tick_params(labelright=False)
b.hist2d(FT+"0",BT+"0",bins=(400,400),range=(lim,lim),log=True)
f.subplots_adjust()

d = .015  # how big to make the diagonal lines in axes coordinates
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
plt.ylim((2620,3450))
f.text(0.5, 0.04, 'TDC channel (front)', ha='center')

a=ax.get_xticks().tolist()
a = [float(b)/1e3 for b in a]
ax.set_xticklabels(a)

plt.savefig("tdc/tdc_raw.pdf",dpi=600)
plt.savefig("tdc/tdc_raw.pgf",dpi=600)
#
# fig tdc_raw
#




exps = ["BT_"+str(i) for i in range(3)]
center = r.asNumpy(exps)
center = [center[e] for e in exps]
center = np.concatenate(center)
center = center[center>1000]
center = np.median(center)
print("center",center)

def get_data(id,strip,side):
    l = []
    for i in range(3):
        a = r.filter("id["+str(i)+"]=="+str(id))
        if "F" in side:
            a = a.filter("FI["+str(i)+"]=="+str(strip))
        else:
            a = a.filter("BI["+str(i)+"]=="+str(strip))
        exp = side+str(i)
        data = a.asNumpy([exp])
        l.append(data[exp])
    return np.concatenate(l)

@anycache(cachedir='merged/cache')
def get_offset(id,strip,side,plot=False):
    result = get_data(id,strip,side)
    if result.sum() == 0: return 0

    lim = (6.5e7,6.6e7)
    lim = (center-5e5,center+5e5)
    count, bins = np.histogram(result,bins=200,range=lim)
    bins = (bins[0:-1]+bins[1:])/2

    def gauss2(p,x):
        a = p[0:2]
        b = p[2:4]
        c = p[4:6]
        return a[0]*np.exp(-np.power(x-b[0],2)/c[0]) + a[1]*np.exp(-np.power(x-b[1],2)/c[1])

    median = np.median(result)
    
    index = np.searchsorted(bins,median)-1
    print(index,bins[index])
    
    def ismax(i):
        return all([count[i]>=count[i+d] for d in range(-8,9)]) and count[i]>count.sum()/600

    lower = 20
    while not ismax(lower):
        lower += 1
    print("lower",lower,bins[lower])
    upper = 180
    while not ismax(upper):
        upper -= 1
    print("upper",upper,bins[upper])

    b_max = bins[upper]
    a_max = count[upper]
    b_min = bins[lower]
    a_min = count[lower]
    p2 = [a_min,a_max,b_min,b_max,25e7,25e7]

    def get_cost_fun(x,y):
        return lambda p: np.power(gauss2(p,x)-y,2).sum()

    f = get_cost_fun(bins,count)
    
    limit = [(a_min*0.7,a_min*1.5),(a_max*0.7,a_max*1.5)]+[(b_min-5e4,b_min+5e4),(b_max-5e4,b_max+5e4)]+[(None,None)]*2
    m = Minuit.from_array_func(f,p2,limit=limit)
    m.migrad()

    if plot:
        x = np.arange(*lim,1000)
        plt.hist(result,bins=200,range=lim)
        plt.plot(x,gauss2(p2,x))
        p = m.np_values()
        plt.plot(x,gauss2(p,x))

    return m.np_values()[2:4].sum()/2

def tovector(x):
    x = [float(xi) for xi in x]
    return 'std::vector<float>{'+str(x).replace('[','').replace(']','')+'}'

L = (1.3e7,1.5e7)
strips = (0,33)
plt.figure(figsize=(5,8))

## DET 0
plt.clf()
get_offset(0,9,BT,plot=True)
plt.savefig("tdc/plot.pdf")
if True:
    offsets = [get_offset(0,i,BT) for i in range(1,25)]
    r = r.define("id0B",tovector(offsets))
    for j in range(3):
        i = str(j)
        r = r.define(BT+"_"+i,"id["+i+"]==0 ? "+BT+i+"-id0B[BI["+i+"]-1] : "+BT+i)
    BT = BT+"_"

## DET 1
if True:
    offsets = [get_offset(1,i,BT) for i in range(1,25)]
    r = r.define("id1B",tovector(offsets))
    for j in range(3):
        i = str(j)
        r = r.define(BT+"_"+i,"id["+i+"]==1 ? "+BT+i+"-id1B[BI["+i+"]-1] : "+BT+i)
    BT = BT+"_"

## DET 2
if True:
    offsets = [get_offset(2,i,BT) for i in range(1,17)]
    r = r.define("id2B",tovector(offsets))
    offsets = [get_offset(2,i,FT) for i in range(1,17)]
    r = r.define("id2F",tovector(offsets))
    for j in range(3):
        i = str(j)
        r = r.define(BT+"_"+i,"id["+i+"]==2 ? "+BT+i+"-id2B[BI["+i+"]-1] : "+BT+i)
        r = r.define(FT+"_"+i,"id["+i+"]==2 ? "+FT+i+"-id2F[FI["+i+"]-1] : "+FT+i)
    BT = BT+"_"
    FT = FT+"_"

    offsets = []
    for i in range(1,17):
        b = r.filter("id[1]==2").define("delta","("+BT+"1 - " + FT+"1)")
        arr = b.filter("BI[1]=="+str(i)).filter("abs(delta)<20e3").asNumpy(["delta"])["delta"]
        offsets = offsets + [arr.mean()]
    r = r.define("id2B2",tovector(offsets))
    for j in range(3):
        i = str(j)
        r = r.define(BT+"_"+i,"id["+i+"]==2 ? "+BT+i+"-id2B2[BI["+i+"]-1] : "+BT+i)
    BT = BT+"_"
    
    offsets = []
    for i in range(1,17):
        b = r.filter("id[1]==2").define("delta","("+BT+"1 - " + FT+"1)")
        arr = b.filter("FI[1]=="+str(i)).filter("abs(delta)<10e3").asNumpy(["delta"])["delta"].mean()
        offsets = offsets + [0 if np.isnan(arr) else -arr]
    print(offsets)
    r = r.define("id2F2",tovector(offsets))
    for j in range(3):
        i = str(j)
        r = r.define(FT+"_"+i,"id["+i+"]==2 ? "+FT+i+"-id2F2[FI["+i+"]-1] : "+FT+i)
    FT = FT+"_"

## DET 3
if True:
    offsets = [get_offset(3,i,BT) for i in range(1,17)]
    r = r.define("id3B",tovector(offsets))
    offsets = [get_offset(3,i,FT) for i in range(1,17)]
    r = r.define("id3F",tovector(offsets))
    for j in range(3):
        i = str(j)
        r = r.define(BT+"_"+i,"id["+i+"]==3 ? "+BT+i+"-id3B[BI["+i+"]-1] : "+BT+i)
        r = r.define(FT+"_"+i,"id["+i+"]==3 ? "+FT+i+"-id3F[FI["+i+"]-1] : "+FT+i)
    BT = BT+"_"
    FT = FT+"_"

## Correction
for j in range(3):
    i = str(j)
    r = r.define("thetaa"+i,"abs(thetaLab["+i+"]-90)")
    r = r.define("dist"+i,"thetaa"+i+"<40 ? \
            sqrt(pow(4,2) + pow(4*tan(thetaa"+i+"/180*3.141592),2)) : \
            sqrt(pow(3.6,2) + pow(3.6*tan((90-thetaa"+i+")/180*3.141592),2))")
    r = r.define("correction"+i,"1e12*dist"+i+"/100/sqrt(2*pow((2.99792458*1e8),2)*eLab["+i+"]*1e-6/(3.727379378))")
    r = r.define(BT+"_"+i,BT+i+"-correction"+i)
    r = r.define(FT+"_"+i,FT+i+"-correction"+i)
BT = BT+"_"
FT = FT+"_"

L = (-5e5,5e5)
for i in range(4):
    j = 2
    for j in range(3):
        plt.subplot(4,3,i*3+j+1)
        b = r.filter("id["+str(j)+"]=="+str(i))
        b.hist2d([FT+str(j)],["FI["+str(j)+"]"],bins=(33,500),range=(L,strips))


r = r.define("dt","max("+BT+"0,max("+BT+"1,"+BT+"2))-min("+BT+"0,min("+BT+"1,"+BT+"2))")
r = r.define("dt_o","max(BT[0],max(BT[1],BT[2]))-min(BT[0],min(BT[1],BT[2]))")

plt.savefig("tdc_compare.pdf",dpi=600)

if True:
    plt.clf()
    plt.figure(figsize=(width,width/2))#/1.61803398875))
    f, (ax, ax2) = plt.subplots(1, 2,sharey=True,figsize=(width,width/2))#/1.61803398875))
    
    i = str(1)
    
    plt.sca(ax)
    b = r.filter("id["+i+"]==2")\
            .define("delta","(FT_"+i+" - BT_"+i+")/1e3")\
            .filter("abs(delta)<1e4")
    b.hist2d("delta","BI["+i+"]",bins=(16,200),range=((90,130),(1,16)))
    plt.xlabel("Front/Back difference [ns]")
    plt.ylabel("Back strip id")
    #a=ax.get_xticks().tolist()
    #print(a)
    #a = [float(b)-100 for b in a]
    #print(a)
    #ax.set_xticklabels(a)
    #b.hist("delta",bins=500)
    ax.set_xticks([90,100,110,120,130])
    ax.set_xticklabels([-20,-10,0,10,20])

    
    plt.sca(ax2)
    b = r.filter("id["+i+"]==2")\
            .define("delta","("+FT+i+" - "+BT+i+")/1e3")\
            .filter("abs(delta)<1e4")
    b.hist2d("delta","BI["+i+"]",bins=(16,200),range=((-20,20),(1,16)))
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.xlabel("Front/Back difference [ns]")
    #plt.ylabel("Back strip id")
    


    plt.savefig("tdc/front_back.pdf",dpi=600)
    plt.savefig("tdc/front_back.pgf",dpi=600)

if True:
    plt.clf()
    plt.figure(figsize=(width,width/2))#1.61803398875))
    #f, (ax, ax2) = plt.subplots(1, 2,figsize=(width,width/1.61803398875))
    #plt.subplot(1,2,1)
    #plt.sca(ax2)

    b = r
    b = b.define("last","max("+BT+"0,max("+BT+"1,"+BT+"2))/1e3")
    b = b.define("first","min("+BT+"0,min("+BT+"1,"+BT+"2))/1e3")
    lim = (0,3e2)
    b.hist2d("first","last",bins=(200,200),range=((0,100),(0,200)))
    plt.plot([0,30e5],[10,30e5+10],'r')
    plt.xlim(lim)
    plt.ylim(lim)
    plt.xlim((0,1e2))
    plt.ylim((0,1.5e2))
    plt.xlabel("First arrival [ns]")

    #plt.subplot(1,2,2)
    #plt.sca(ax)
    #b = r_back
    #b = b.define("last_","max(BT_0,max(BT_1,BT_2))")
    #b = b.define("first_","min(BT_0,min(BT_1,BT_2))")
    #data = b.asNumpy(["last_","first_"])
    #offset = data["first_"][data["first_"]>10e5].min()
    #offset2 = data["last_"][data["last_"]>10e5].min()
    #b = b.define("last","(last_-"+str(offset2)+")/1e3-0.7e3")
    #b = b.define("first","(first_-"+str(offset)+")/1e3-0.72e3")
    #lim = (0,300)
    #lim2 = (0,150)
    #lim = (0,300)
    #lim2 = (0,450)
    #b.hist2d("first","last",bins=(200,100),range=(lim,lim2),log=True)
    #plt.xlim(lim)
    #plt.ylim(lim2)
    #b.hist("first",bins=200,range=(0,25e2))
    #b.hist("last",bins=200,range=(0,25e2))
    plt.ylabel("Last arrival [ns]")
    #plt.xlabel("First arrival [ns]")


    plt.savefig("tdc/comparison.pdf",dpi=600)
    plt.savefig("tdc/comparison.pgf",dpi=600)


plt.clf()

r = r.define("esum","eCM[0]+eCM[1]+eCM[2]")
r = r.define("X","((eCM[0]+2*eCM[1])/esum-1)/sqrt(3)")
r = r.define("Y","eCM[0]/esum-1.0/3.0")
r = r.define("rho_","sqrt(pow(X,2)+pow(Y,2))*3")
#r = r.filter("!(dt_o>50e3 && dt_o<100e3 && dt>30e3)")

plt.figure(figsize=(15,7))
v = ["dt","dt_o"]
title = ["Calibrated","Uncalibrated"]
for i in range(2):
    plt.subplot(1,2,i+1)
    plt.title(title[i])
    plt.xlabel("$\\Delta t\\ [\\mathrm{ps}]$")
    plt.ylabel("$E_\\mathrm{Lab}\\ [\\mathrm{keV}]$")
    r.hist2d(v[i],"rho_",bins=(500,500),range=((0,300e3),(0,1.0)))
    rho = np.linspace(0,1.0,100)
    plt.plot(np.power(rho,5)*5e3+10e3,rho,'r')

plt.savefig("tdcfit.pdf",dpi=600)

#r = r.filter("dt<100e3")
filt = "(dt<10e3 || dt>1e6)"
filt = "(dt<(10e3+pow(rho_,5)*5e3) || dt>1e6)"
#filt = "(dt<(10e3+pow(rho_,5)*5e3)) "
b = r
#r = r.filter("(dt<(20e3+3e5*pow(rho,2)) || dt>5000e3)")
r = r.filter(filt)
# eff
c = r
c.hist(["eLab[0]","eLab[1]","eLab[2]"],bins=200,range=(0,800),export=("tdc/eff/testname_.npz","count"))
for i in range(3):
    c2 = c.filter("BT["+str(i)+"]==0")
    c2.hist("eLab["+str(i)+"]",bins=200,range=(0,800),export=("tdc/eff/testname_"+str(i)+".npz","count"))

# eff
r.snapshot("a","output/data_tdc_cut.root",["dt","eCM","eLab","mul","pT","deltaE","theta*","phi*","exC12","BI","FI","id",BT+"0",BT+"1",BT+"2","BT_0","BT_1","BT_2","rho_","dist*","dist0","dist1","dist2","correction0","correction1","correction2"])


