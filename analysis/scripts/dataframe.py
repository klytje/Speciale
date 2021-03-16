
from turbo import *
import matplotlib.pyplot as plt
import numpy as np
from ROOT import RDataFrame, TChain
import ROOT
import os
from hashlib import md5
ROOT.ROOT.EnableImplicitMT()
from matplotlib.colors import LogNorm

class dataframe:
    data = None
    def __init__(self,chain,data=None):
        if data is None:
            self.data = RDataFrame(chain)
        else:
            self.data = data

    def filter(self,s):
        return dataframe(None,data=self.data.Filter(s))

    def define(self,s1,s2):
        return dataframe(None,data=self.data.Define(s1,s2))

    def asNumpy(self,names):
        def transform(name):
            tname = 'x'+md5(name.encode('utf-8')).hexdigest()
            cols = self.data.GetColumnNames()
            if name in cols:
                return name
            elif tname not in cols:
                self.data = self.data.Define(tname,name)
            return tname
        for n in names: transform(n)
        res = self.data.AsNumpy([transform(n) for n in names])
        for n in names:
            if transform(n) in res.keys():
                res[n] = res.pop(transform(n))
        return res

    def export(self,name,data,fname):
        if "npz" not in fname: fname = fname + ".npz"
        d = dict()
        if os.path.exists(fname):
            f = np.load(fname)
            for n in f:
                d[n] = f[n]
        d[name]=data
        np.savez_compressed(fname,**d)

    def hist(self,x,bins=None,range=None,label=None,density=False,color=None,export=None,weights=[]):
        if type(x) is str:
            x = [x]
        if type(weights) is str:
            weights = [weights]
        data = self.asNumpy(x+weights)
        d = np.concatenate([data[i] for i in x])
        if len(weights)>0:
            w = np.concatenate([data[i] for i in weights])
        else:
            w = None

        count,bins = np.histogram(d,bins=bins,range=range,density=density,weights=w)
        plt.hist(bins[:-1],histtype='step',bins=bins,weights=count,label=label,range=range\
                ,density=density,color=color)

        if export is not None:
            self.export(export[1],count,export[0])
            self.export("bins",bins,export[0])

    def hist2d(self,x,y,bins=10,range=None,export=None,weights=[],log=False,vmax=None):
        if type(x) is str: x = [x]
        if type(y) is str: y = [y]
        if type(weights) is str: weights = [weights]
        d = self.asNumpy(x+y+weights)
        dx = np.concatenate([d[i] for i in x])
        dy = np.concatenate([d[i] for i in y])
        if len(weights)>0:
            #w = self.asNumpy(weights)
            w = np.concatenate([d[i] for i in weights])
        else:
            w = None
        
        #print(w)
        H, xedges, yedges = np.histogram2d(dy,dx,range=(range[1],range[0]),bins=bins,weights=w)

        if export is not None:
            self.export(export[1],H.flatten(),export[0])
            self.export("bins",(xedges,yedges),export[0])

        mpl_data = RGBToPyCmap(turbo_colormap_data)
        plt.register_cmap(name='turbo', data=mpl_data, lut=turbo_colormap_data.shape[0])
        H[H == 0] = np.nan
        if not log:
            c = plt.imshow(H[::-1,:], extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),vmin=0,aspect='auto',vmax=vmax,cmap='turbo')
        else:
            c = plt.imshow(H[::-1,:], extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]),vmin=1,aspect='auto',norm=LogNorm(vmin=1,vmax=vmax),cmap='turbo')
        plt.colorbar(c)
    
    def getColumnNames(self):
        return self.data.GetColumnNames()

    def snapshot(self,treename,filename,columns=[]):
        if columns:
            columns = "|".join(columns)
        else:
            columns = ""
        return self.data.Snapshot(treename,filename,columns)
