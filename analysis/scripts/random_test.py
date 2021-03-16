import numpy as np
from matplotlib import pyplot as plt
import sys
from mcfit import mcfitter
from datatools import *

folder = "interference"
folder2 = "interference_j1"

data,bins = load_data("data")
balamuth,w = get_weighted(folder)
balamuth_j1,w_j1 = get_weighted(folder2)



def getrandom(folder,r,delta):
    rand, bins = load_data("random")
    _, eCM_bins = load_data("random","eCM_predet")
    eCM_bins = (eCM_bins[0:-1]+eCM_bins[1:])/2

    H = np.reshape(rand,(100,100))

    yedges = bins[0]
    xedges = bins[1]

    X = (xedges[0:100]+xedges[1:101])/2
    Y = (yedges[0:100]+yedges[1:101])/2
    X,Y = np.meshgrid(X,Y)

    K = 9100.0

    E1 = K/3*(3*Y+1)
    E2 = 1.0/6*(3*np.sqrt(3)*K*X-3*K*Y+2*K)
    E3 = 1.0/6*(-3*np.sqrt(3)*K*X-3*K*Y+2*K)


    def get(E,k,delta):
        balamuth_E, w_E = get_weighted(folder,"eCM_predet")
        s = np.interp(E,eCM_bins,balamuth_E*w_E(k,delta))
        s = s/s.sum()
        return s

    def weight(k,delta):
        s = 1
        s *= get(E1,k,delta)
        s *= get(E2,k,delta)
        s *= get(E3,k,delta)
        return s
    
    return H,weight

H, weight_ = getrandom("interference",0.296,0.6757)

H[H==0] = np.nan
#c = plt.imshow(H[::-1,:],extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),aspect='auto',vmin=0)
yedges = bins[0]
xedges = bins[1]
c = plt.imshow((weight_(0.296,0.6757)*H)[::-1,:],extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]),aspect='auto',vmin=0)

plt.savefig("random.pdf")
