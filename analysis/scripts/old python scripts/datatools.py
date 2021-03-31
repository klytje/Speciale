import numpy as np
from matplotlib import pyplot as plt

def fix(a,balamuth):
    a = a/balamuth
    a[np.isnan(a)] = 0
    return a

def load_data(name,figure="dalitz_forfit",path=""):
    if path == "":
        bins = np.load("figures/"+name+"/"+figure+".npz")
    else:
        bins = np.load(path+"/"+figure+".npz")
    data = bins["count"]
    bins = bins["bins"]
    return data,bins

def get_weighted(folder,figure="dalitz_forfit",path=""):
    if path == "":
        data = np.load("figures/"+folder+"/"+figure+".npz")
    else:
        data = np.load(path+'/'+folder+'/'+figure+'.npz')
    balamuth = data["count"]
    f1 = fix(data["f1"],balamuth)
    f2 = fix(data["f2"],balamuth)
    re = fix(data["re"],balamuth)
    im = fix(data["im"],balamuth)
    w = lambda k,d: k*f1+(1-k)*f2\
            +2*np.sqrt(k*(1-k))*(re*np.cos(2*np.pi*d)+im*np.sin(2*np.pi*d))
    return balamuth.astype(float),w

def plot(mcf,name):
    print("chi2/ndf: ",str(int(mcf.chisquare()))+"/"+str(mcf.ndf()))
    print("pval: ",str(mcf.pvalue()))
    def subplot(H,title):
        norm = np.max(np.abs(H))
        xedges = mcf.bins[0]
        yedges = mcf.bins[1]
        #H[H==0] = np.nan
        c = plt.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',cmap='seismic',vmin=-norm,vmax=norm)
        plt.gca().set_aspect("equal")
        plt.xlabel("$X$")
        plt.ylabel("$Y$")
        plt.title(title)
        plt.colorbar(c)
    
    plt.clf()
    plt.figure(figsize=(15,7))

    plt.subplot(1,2,1)
    H = np.reshape(mcf.data-mcf.prediction(),(100,100))
    subplot(H,"$\\mathrm{Residuals (data-fit)}$")

    plt.subplot(1,2,2)
    H = np.reshape(np.sign(mcf.data-mcf.prediction())*mcf.chisquare_binwise(),(100,100))
    subplot(H,"$\\chi^2 \\cdot \\mathrm{sign(data-fit)}$")
    plt.savefig(name+".pdf")

